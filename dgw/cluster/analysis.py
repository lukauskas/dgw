from collections import defaultdict
from logging import debug
from math import floor
import scipy.cluster.hierarchy as hierarchy
import pandas as pd
import numpy as np
from ..data.containers import AlignmentsData
from ..dtw.distance import dtw_std, dtw_path_is_reversed, warping_conservation_vector
from ..dtw import transformations, dtw_projection, no_nans_len
from ..dtw.parallel import parallel_dtw_paths
import gzip
import scipy.stats

from scipy.cluster._hierarchy import get_max_dist_for_each_cluster

def compute_paths(data, dtw_nodes_list, n, n_processes=None, *dtw_args, **dtw_kwargs):

    non_leaf_nodes = dtw_nodes_list[n:]
    paths = parallel_dtw_paths(data, non_leaf_nodes, n_processes=n_processes, *dtw_args, **dtw_kwargs)
    return paths

def _to_dtw_tree(linkage, hierarchical_clustering_object, prototypes, prototyping_function='mean'):
    """
    Converts a hierarchical clustering linkage matrix `linkage` to hierarchy of `DTWClusterNode`s.
    This is a modification of `scipy.cluster.hierarchy.to_tree` function and the code is mostly taken from it.

    :param linkage: linkage matrix to convert to the DTW Tree
    :param hierarchical_clustering_object: hierarchical clustering object to work with
    :param prototyping_function: "reduce" function for prototype calculation, or "mean" to simply use data mean
    """

    # Validation
    linkage = np.asarray(linkage, order='c')
    hierarchy.is_valid_linkage(linkage, throw=True, name='Z')

    data = hierarchical_clustering_object.data
    labels = data.items
    values = data.ix

    n = linkage.shape[0] + 1

    # Create a list full of None's to store the node objects
    d = [None] * (n * 2 - 1)

    # Create the nodes corresponding to the n original objects.
    for i in xrange(0, n):
        index = labels[i]
        d[i] = DTWClusterNode(id=index, hierarchical_clustering_object=hierarchical_clustering_object,
                              prototype=values[index])

    nd = None

    for i in xrange(0, n - 1):
        fi = int(linkage[i, 0])
        fj = int(linkage[i, 1])

        assert(fi <= i + n)
        assert(fj <= i + n)

        id = i + n
        left = d[fi]
        right = d[fj]
        dist = linkage[i, 2]

        if prototypes:
            prototype = prototypes[id]

            nd = DTWClusterNode(id=id, hierarchical_clustering_object=hierarchical_clustering_object,
                                prototype=prototype,
                                left=left, right=right,
                                dist=linkage[i, 2])

        elif callable(prototyping_function):
            prototype = prototyping_function(left.prototype.values, right.prototype.values, left.count, right.count)

            nd = DTWClusterNode(id=id, hierarchical_clustering_object=hierarchical_clustering_object,
                                prototype=prototype,
                                left=left, right=right,
                                dist=linkage[i, 2])

        elif prototyping_function == 'mean':
            nd = DTWClusterNode(id=id, hierarchical_clustering_object=hierarchical_clustering_object,
                                prototype=None,
                                left=left, right=right,
                                dist=linkage[i, 2])

            # A bit hacky, but does job. Doing this as to get to use nd.data
            nd._prototype = nd.data.mean()

        assert(linkage[i, 3] == nd.count)
        d[n + i] = nd

    return nd, d

def add_path_data(dtw_nodes, n, paths):
    """
    Adds precomputed path data to dtw_nodes
    :param data:
    :param dtw_nodes:
    :param n:
    :param paths:
    :return:
    """
    # Loop through non-leaf nodes
    for node in dtw_nodes[n:]:
        node.warping_paths = paths[node.id]

# Inheritting from object explicitly as hierarchy.ClusterNode is not doing this
class DTWClusterNode(object, hierarchy.ClusterNode):
    _prototype = None
    _index = None
    _hierarchical_clustering_object = None
    _projected_data = None
    _warping_paths = None
    _tracked_points = None
    _tracked_points_histogram = None
    _points_of_interest_histogram = None
    _warping_conservation_data = None

    def __init__(self, hierarchical_clustering_object, id, prototype, left=None, right=None, dist=0, count=1):
        hierarchy.ClusterNode.__init__(self, id, left=left, right=right, dist=dist, count=count)
        if not isinstance(prototype, pd.DataFrame):
            prototype = pd.DataFrame(prototype, columns=hierarchical_clustering_object.data.dataset_axis)
        self._prototype = prototype
        self._hierarchical_clustering_object = hierarchical_clustering_object
        self._index = pd.Index(self.__get_item_ids())

        # Assume no points of interest. TODO: Add  way to specify those
        self._points_of_interest = {}

    def __get_item_ids(self):
        """
        Gets the ids of all the leaf nodes that are in the tree below
        :return:
        """
        if self.is_leaf():
            return pd.Index([self.id])
        else:
            return self.get_left().index | self.get_right().index

    @property
    def prototype(self):
        return self._prototype

    @property
    def index(self):
        return self._index

    def reindex(self, new_index):
        # Make sure they are the same elements just in different order
        assert(len(new_index & self._index) == len(self._index) == len(new_index))
        self._index = new_index
        return self._index

    @property
    def data(self):
        # Don't store the data inside
        return self._hierarchical_clustering_object.data.ix[self.index]

    @property
    def projected_data(self):
        # It is infeasible to calculate projected data for all the nodes beforehand
        self.ensure_projections_are_calculated()
        return self._projected_data




    def ensure_projections_are_calculated(self):
        if not self._projected_data:
            self._projected_data = self.__project_items_onto_prototype()

    def _calculate_histogram(self, points_of_interest, number_of_bins, lengths=None):

        histogram = defaultdict(lambda: np.zeros(number_of_bins))
        for ix, poi in points_of_interest.iteritems():
            if lengths is not None:
                scaling_ratio = float(number_of_bins) / lengths[ix]
            else:
                scaling_ratio = 1

            for poi_name, points in poi.iteritems():
                current_histogram = histogram[poi_name]
                for point in set(points):
                    min_rescaled = int(floor(point * scaling_ratio))
                    max_rescaled = int(floor((point + 1) * scaling_ratio))
                    for rescaled_point in xrange(min_rescaled, max_rescaled):
                        assert rescaled_point <= number_of_bins
                        current_histogram[rescaled_point] += 1

        return pd.DataFrame(histogram)

    @property
    def tracked_points_histogram(self):
        if self._tracked_points_histogram is None:
            self.ensure_points_of_interest_are_tracked_down()
            self._tracked_points_histogram = self._calculate_histogram(self._tracked_points, len(self.prototype))
        return self._tracked_points_histogram

    @property
    def points_of_interest_histogram(self):
        if self._points_of_interest_histogram is None:
            self._points_of_interest_histogram = self._calculate_histogram(self.points_of_interest,
                                                                           max(self.data.lengths),
                                                                           self.data.lengths)
        return self._points_of_interest_histogram

    def ensure_points_of_interest_are_tracked_down(self):
        if self._tracked_points is None:
            self._tracked_points = self.__track_points_of_interest()


    @property
    def warping_paths(self):
        if not self._warping_paths:
            raise Exception('Warping paths not computed')

        return self._warping_paths

    @warping_paths.setter
    def warping_paths(self, values):
        self._warping_paths = values

    @property
    def regions(self):
        """
        :rtype: Regions
        """
        if self._hierarchical_clustering_object.regions is None:
            return None
        else:
            return self._hierarchical_clustering_object.regions.ix[self.index]

    def __compute_dtw_warping_paths(self):
        data = self.data
        print(data)
        dtw_function = self._hierarchical_clustering_object.dtw_function
        prototype = self.prototype

        paths = {}
        for ix in data.items:
            _, _, path = dtw_function(data.ix[ix], prototype, dist_only=False)
            # Reduce the bit size of the path arrays to 16 bits
            # DTW would be too slow to use anyway if we had more than 2**16-1 items in it
            # Feel free to update this if it is not the case.
            #path = (np.asarray(path[0], dtype=np.int16), np.asarray(path[1], dtype=np.int16))

            paths[ix] = path

        return paths

    def save_as_encode_track(self, filename, track_name=None, track_description=None):
        if self.regions is None:
            raise Exception('Cannot save {0!r} as region locations are not specified'.format(self))

        if track_name:
            track_kwargs = {'name': track_name}
            if track_description:
                track_kwargs['description'] = track_description
        else:
            track_kwargs = {}

        regions = self.regions

        regions = regions.infer_strand_from_whether_the_region_was_reversed_or_not(self.reversal_dictionary)
        regions.to_bed(filename, **track_kwargs)

    def save_prototype_to_text(self, filename):

        prototype = self.prototype
        f = open(filename, 'w')
        try:

            f.write('#{0}'.format('bin'))
            for col in prototype.columns:
                f.write('\t{0}'.format(col))
            f.write('\n')

            for bin, row in prototype.iterrows():
                f.write('{0}'.format(bin))
                for col in prototype.columns:
                    f.write('\t{0}'.format(row[col]))
                f.write('\n')
        finally:
            f.close()

    def save_conservation_coefficient_as_text(self, filename):

        conservation_vector = self.warping_conservation_vector()

        f = open(filename, 'w')
        try:
            f.write('#{0}\t{1}\t{2}\n'.format('start_bin', 'end_bin', 'avg_conservation'))

            for i in xrange(len(conservation_vector)):
                f.write('{0}\t{1}\t{2}\n'.format(i, i+1, conservation_vector[i]))
        finally:
            f.close()


    def save_as_list_of_indices(self, filename):
            index = self.data.items
            f = open(filename, 'w')
            try:
                for ix in index:
                    f.write('{0}\n'.format(ix))
            finally:
                f.close()

    def save_poi_histograms_to_file(self, basename):
        points_of_interest = self.points_of_interest
        if not points_of_interest:
            return

        raw_filename = basename + '-raw.csv'
        warped_filename = basename + '-warped.csv'

        self.points_of_interest_histogram.to_csv(raw_filename)
        self.tracked_points_histogram.to_csv(warped_filename)

    def poi_entropies(self):
        if not self.points_of_interest:
            return None

        untracked_hist = self.points_of_interest_histogram
        tracked_hist = self.tracked_points_histogram

        untracked_entropy = pd.Series(scipy.stats.entropy(untracked_hist),
                                      index=untracked_hist.columns, name='raw')
        tracked_entropy = pd.Series(scipy.stats.entropy(tracked_hist),
                                    index=tracked_hist.columns, name='warped')

        df = pd.DataFrame([untracked_entropy, tracked_entropy]).T
        df.index.name = 'poi_file'

        df['diff'] = df['raw'] - df['warped']
        df['rel_diff'] = df['diff'] / df['raw']

        return df

    def save_pois_to_file(self, filename):
        points_of_interest = self.points_of_interest
        if not points_of_interest:
            return

        with gzip.GzipFile(filename, 'w') as f:
            f.write('#region\tpoi_file\tbins\tprototype_bins\n')

            points_of_interest = self.points_of_interest
            warped_points_of_interest = self.tracked_points_of_interest

            for region, poi_data in points_of_interest.iteritems():
                warped_poi_data = warped_points_of_interest[region]

                for poi_filename, pois in poi_data.iteritems():
                    warped_pois = warped_poi_data[poi_filename]

                    str_pois = ';'.join(map(str, pois))
                    str_warped_pois = ';'.join(map(str, warped_pois))

                    f.write('{}\t{}\t{}\t{}\n'.format(region, poi_filename, str_pois, str_warped_pois))


    def save_warpings_to_file(self, filename):
        data = self.data
        index = data.items
        regions = self.regions
        if not self.is_leaf():
            warping_paths = self.warping_paths
        else:
            # if we're at a leaf node, there are no warping paths
            # so we just create sample warping paths that just map each point in data to itself

            assert len(data.items) == 1  # assumption for sanity
            first_ix = data.items[0]
            data_len = len(data.ix[first_ix].dropna())
            warping_paths = {first_ix: np.vstack([np.arange(data_len),
                                                  np.arange(data_len)])}

        if regions:
            bin_intervals = regions.ix[index].bins_to_intervals(data.resolution)
            chromosomes = regions.chromosome


            f = gzip.open(filename, 'w')
            f.write('#{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format('index', 'bin', 'chromosome', 'start', 'end', 'prototype_bin'))
            try:
                for ix, chromosome in chromosomes.iteritems():
                    path = warping_paths[ix]
                    bi = bin_intervals[ix]

                    for i in xrange(len(path[0])):
                        p_a = path[0][i]
                        p_b = path[1][i]

                        current_bin = bi[p_a]

                        f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(ix, p_a, chromosome,
                                                                        current_bin[0], current_bin[1],
                                                                        p_b))
            finally:
                f.close()
        else:
            f = gzip.open(filename, 'w')
            try:
                f.write('#{0}\t{1}\t{2}\n'.format('index', 'bin', 'prototype_bin'))

                for ix, path in warping_paths.items():
                    for i in xrange(len(path[0])):
                        p_a = path[0][i]
                        p_b = path[1][i]
                        f.write('{0}\t{1}\t{2}\n'.format(ix, p_a, p_b))
            finally:
                f.close()


    def save_warping_conservation_data_to_file(self, filename):
        data = self.data
        index = data.items
        conservation_data = self.warping_conservation_data

        f = gzip.open(filename, 'w')
        f.write('#{0}\t{1}\t{2}\t{3}\n'.format('index', 'start_bin', 'end_bin', 'conservation_coefficient'))
        try:
            for ix in index:
                conservation_vector = conservation_data.ix[ix]
                start_i = None
                current_value = None
                for i in xrange(len(conservation_vector)):
                    value = conservation_vector[i]
                    if value < 1e-6:
                        if start_i is not None: # 1e-6 to deal with floating-point issues
                            # If were at zero, some conserved region finishes here
                            f.write('{0}\t{1}\t{2}\t{3}\n'.format(ix, start_i, i, current_value))
                            start_i = None
                            current_value = None
                        else:
                            pass # Nothing to do here
                    elif start_i is None:
                        start_i = i
                        current_value = value
                    else:
                        try:
                            assert(abs(value - current_value) < 1e-6)
                        except AssertionError:
                            debug(value, current_value)
                            raise
        finally:
            f.close()


    def __project_items_onto_prototype(self):

        if self.is_leaf():
            # there is no projection, really, return itself
            return self.data
        else:
            data = self.data
            prototype = self.prototype

            warping_paths = self.warping_paths

            columns = data.dataset_axis

            projections = {}
            for ix in data.items:
                item = data.ix[ix]

                projection = dtw_projection(item, prototype, path=warping_paths[ix])
                df = pd.DataFrame(projection, index=range(len(prototype)), columns=columns)
                projections[ix] = df
            panel = pd.Panel(projections)
            panel = panel.ix[self.data.items]
            ad = AlignmentsData(panel, self.data.resolution)
            return ad

    def _compute_warping_conservation_data(self):
        data = self.data
        if self.is_leaf():
            conservation_data = [np.ones(no_nans_len(self.prototype)-1)]
        else:
            warping_paths = self.warping_paths
            conservation_data = []
            for ix in data.items:
                path = warping_paths[ix]
                conservation_data.append(warping_conservation_vector(path))

        return pd.DataFrame(conservation_data, index=self.index)

    @property
    def warping_conservation_data(self):
        if self._warping_conservation_data is None:
            self._warping_conservation_data = self._compute_warping_conservation_data()
        return self._warping_conservation_data

    def warping_conservation_vector(self):
        conservation_vector = self.warping_conservation_data.mean()
        return conservation_vector

    def __track_points_of_interest(self):

        if self.is_leaf():
            # Nothing to track
            return self.points_of_interest
        else:
            points_of_interest = self.points_of_interest
            warping_paths = self.warping_paths
            tracked_points = defaultdict(lambda: {})
            for ix, pois in points_of_interest.iteritems():
                for j, poi in pois.iteritems():
                    mapped_points = transformations.points_mapped_to(poi, warping_paths[ix])
                    tracked_points[ix][j] = mapped_points

            return tracked_points

    @property
    def points_of_interest(self):
        poi = self.data.points_of_interest
        return poi

    @property
    def reversal_dictionary(self):
        if self.is_leaf():
            # If we are leaf node, there is no DTW projection onto prototype
            # and therefore strand cannot be inferred.
            # Just return None
            return {ix: None for ix in self.index}

        warping_paths = self.warping_paths
        reversal_dict = {}
        for item in self.index:
            reversal_dict[item] = dtw_path_is_reversed(warping_paths[item])

        return reversal_dict

    @property
    def tracked_points_of_interest(self):
        self.ensure_points_of_interest_are_tracked_down()
        tracked_points = self._tracked_points
        return tracked_points

    @property
    def n_items(self):
        """
        Returns number of leaves below the current node
        """
        return len(self.index)

    def __len__(self):
        return self.n_items

    def __repr__(self):
        return "<{0} containing {1} items>".format(self.__class__.__name__, self.n_items)


def _reduce_tree(tree, reduce_func, map_function=lambda node: node.id,
                 is_value=lambda x: not isinstance(x, hierarchy.ClusterNode)):
    """
    Performs reduce operation on the binary tree generated by hierarchical clustering.

    Implementation does not use recursion in order not to overflow the stack

    Equivalent to:
    ```
    def reduce(tree, reduce_function, map_function):
    if tree.is_leaf():
        return map_function(tree)
    else:
        left = reduce(tree.get_left(), reduce_function, map_function)
        right = reduce(tree.get_right(), reduce_function, map_function)
        return reduce_function(tree, left, right)
    ```
    :param tree: root node of the tree
    :type tree: `scipy.cluster.hierarchy.ClusterNode`
    :param reduce_func: function that will be performed on reduction. Takes three parameters: node, left branch, right branch
    :param map_function: function that gets the value of a `ClusterNode` defaults to retrieving the node's id
    :param is_value: function that determines whether a (reduced) item or just another node in cluster
    :return: The result of reduction function
    """

    def _add_children_to_stack(node):
        stack.append(node)  # Append self to stack so we can trace it back later
        stack.append(node.get_left())
        stack.append(node.get_right())

    stack = [tree]

    while stack:
        node = stack.pop()
        if not is_value(node):
            if node.is_leaf():
                try:
                    stack.append(map_function(node))
                except Exception, e:
                    print "Got {0!r} while getting value of {1!r}".format(e, node)
                    raise
            else:
                _add_children_to_stack(node)
        else:
            try:
                next_node = stack.pop()
            except IndexError:
                return node

            if not is_value(next_node):
                stack.append(node)
                stack.append(next_node)
            else:
                parent_node = stack.pop()
                assert(not is_value(parent_node))

                try:
                    reduced_value = reduce_func(parent_node, node, next_node)
                except Exception, e:
                    print "Got {0!r} on reduce of {1!r} {2!r}".format(e, node, next_node)
                    raise

                stack.append(reduced_value)

class HierarchicalClustering(object):
    _condensed_distance_matrix = None
    _data = None
    _linkage_matrix = None
    _distance_threshold = None
    __dtw_args = None
    __dtw_kwargs = None
    _regions = None

    __tree = None
    __tree_nodes_list = None


    def __init__(self, data, regions, linkage_matrix, dtw_function=dtw_std, prototypes=None, prototyping_method='psa'):
        """
        Initialises hierarchical clustering analyser.
        Handles linkage calculation, dendrogram plotting and prototype generation.

        Supports three prototyping methods:

            Prioretised Shape Averaging (PSA) (prototyping_method='psa')
                as described in [#Niennattrakul:2009ep], see also `dgw.transformations.sdtw_averaging`
            Standard average of DTW path (prototyping_method='standard')
                averages the DTW paths for each pair of nodes with the same parent.
                The weights are determined by how many sequences were averaged into the node.
                see `dgw.transformations.dtw_path_averaging`
            Unweighted standard average of DTW path (prototyping_method='standard-unweighted')
                similar to the standard method above, but does not bias sequences higher up the tree
            Mean (prototyping_method='mean')
                simply the mean of the data


        .. [#Niennattrakul:2009ep] Vit Niennattrakul and Chotirat Ann Ratanamahatana "Shape averaging under Time Warping",
           2009 6th International Conference on Electrical Engineering/Electronics, Computer,
           Telecommunications and Information Technology (ECTI-CON)

        :param data: a pd.DataFrame object of the data in clusters
        :type data: AlignmentsData
        :param linkage_matrix: linkage_matrix computed by fastcluster.linkage.
        :param dtw_function: DTW calculation function
        :param prototypes: cluster node prototypes (will be computed if None)
        :param prototyping_method: Averaging method either 'psa', 'standard', 'standard-unweighted', 'mean'
        :return:
        """
        if not isinstance(data, AlignmentsData):
            raise ValueError('Data should be instance of {0}'.format(AlignmentsData.__name__))


        self._data = data
        self._regions = regions

        # small negative distances in linkage matrix are sometimes possible due to rounding errors. Change them to 0
        linkage_matrix[:, 2][linkage_matrix[:, 2] < 0] = 0
        self._linkage_matrix = linkage_matrix

        self.__dtw_function = dtw_function

        tree, tree_nodes = self.__dtw_tree_from_linkage(linkage_matrix, prototypes, prototyping_method)
        self.__tree = tree
        self.__tree_nodes_list = tree_nodes

    def extract_prototypes(self):
        prototypes = {}

        for node in self.tree_nodes_list:
            prototypes[node.id] = node.prototype

        return prototypes

    def __dtw_tree_from_linkage(self, linkage, prototypes,  method):
        """
        Computes a prototyped tree from linkage matrix
        :param linkage: linkage matrix
        :param prototypes: possibly precomputed prototypes
        :param method: prototyping method
        :return:
        """



        if method == 'psa':
            averaging_func = lambda x, y, wx, wy: \
                transformations.sdtw_averaging(x, y, wx, wy, dtw_function=self.dtw_function)
        elif method == 'standard':
            averaging_func = lambda x, y, wx, wy: \
                transformations.dtw_path_averaging(x, y, wx, wy, dtw_function=self.dtw_function)
        elif method == 'standard-unweighted':
            averaging_func = lambda x, y, wx, wy: \
                transformations.dtw_path_averaging(x, y, 1, 1, dtw_function=self.dtw_function)
        elif method == 'mean':
            averaging_func = 'mean' # Not a function really, but the code will deal with it
        else:
            raise ValueError('Incorrect method supplied: '
                             'only \'psa\', \'standard\' or \'standard-unweighted\' supported')

        return _to_dtw_tree(linkage, self, prototypes, averaging_func)

    @property
    def data(self):
        """
        Returns the data that is internal to the object
        :rtype: AlignmentsData
        """
        return self._data

    @property
    def regions(self):
        return self._regions

    @property
    def dtw_function(self):
        return self.__dtw_function

    @property
    def condensed_distance_matrix(self):
        return self._condensed_distance_matrix

    @property
    def linkage(self):
        """
        Returns the linkage matrix of the clustering provided
        :return: linkage matrix.
        """
        return self._linkage_matrix

    def distance_threshold_for_n_clusters(self, n_clusters):
        """
        Returns distance threshold that can be used to forn `n_clusters`.
        :param n_clusters: number of clusters to form
        :return:
        """
        n_clusters = int(n_clusters)

        assert n_clusters > 0, 'Minimum number of clusters is 1, got {}'.format(n_clusters)
        linkage = self.linkage

        n = self.num_obs
        assert n >= n_clusters, 'Specified number of clusters ' \
                                '{} is larger than number of data points {}'.format(n_clusters, n)

        # Special case, otherwise it doesn't work
        if n_clusters == 1:
            return np.inf
        else:
            max_distances = np.empty(self.num_obs, dtype=np.double)
            get_max_dist_for_each_cluster(linkage, max_distances, n)

            threshold = max_distances[-n_clusters]

            return threshold

    @property
    def num_obs(self):
        return self.data.number_of_items

    @property
    def dataset_names(self):
        return self.data.dataset_axis

    def as_tree(self):
        return self.__tree

    @property
    def tree_nodes_list(self):
        return self.__tree_nodes_list

    def dendrogram(self, ax=None, no_labels=True, *args, **kwargs):
        """
        Plots the dendrogram for the hieararchical clustering
        :return:
        """
        from dgw.util.plotting import pyplot as plt

        linkage = self.linkage

        if ax is None:
            ax = plt.gca()

        color_threshold = kwargs.pop('color_threshold', self._distance_threshold)

        ans = hierarchy.dendrogram(linkage, no_labels=no_labels,
                                   above_threshold_color='k',
                                   color_threshold=color_threshold,
                                   *args, **kwargs)

        ax.set_xlabel('Distance')

        return ans



    def _rename_nodes(self, tree):
        """
        Renames the leaf nodes of the cluster tree generated by linkage to match the actual indices of the data.
        :param tree:
        :return:
        """
        stack = [tree]

        while stack:
            node = stack.pop()

            if node.is_leaf():
                node.id = self.data.items[node.id]
            else:
                stack.append(node.get_left())
                stack.append(node.get_right())

    def cut(self, t):
        """
        Cuts the dendrogram at specified threshold t.
        :param t: threshold
           :return:
        """

        self._distance_threshold = t
        root = self.as_tree()

        queue = set([root])

        clusters = set()
        while queue:
            current_node = queue.pop()
            if current_node.dist >= t:
                queue.add(current_node.get_right())
                queue.add(current_node.get_left())
            else:
                clusters.add(current_node)
        return ClusterAssignments(self, clusters, t)

    def cut_and_resort(self, cut_threshold, index):
        """
        Cuts the dendrogram based on the color list returned by `scipy.cluster.hierarchy.dedrogram`
        :param cut_threshold: cut threhsold to cut clusters at
        :param index: index of nodes that is already in order
        :return: clusters
        """
        cluster_assignments = self.cut(cut_threshold)


        already_asigned_indices = pd.Index([])
        for cluster in cluster_assignments:
            cluster_index = cluster.index
            sub_index = pd.Index([i for i in index if i in cluster_index])
            cluster.reindex(sub_index)

            if len(already_asigned_indices & sub_index):
                raise Exception("There is some overlap between cluster cuts. There shouldn't be")
            already_asigned_indices =  (already_asigned_indices | sub_index)

        return cluster_assignments


class ClusterAssignments(object):
    _hierarchical_clustering_object = None
    _cluster_roots = None

    _clusters = None
    _cut_depth = None

    def __init__(self, hierarchical_clustering_object, cluster_roots, cut_depth):
        self._hierarchical_clustering_object = hierarchical_clustering_object
        self._cluster_roots = cluster_roots

        clusters = []

        for cluster_root in self._cluster_roots:
            clusters.append(cluster_root)
        self._clusters = clusters
        self._cut_depth = cut_depth

    @property
    def n(self):
        return len(self._cluster_roots)

    def __len__(self):
        return self.n

    @property
    def dataset_names(self):
        return self.hierarchical_clustering_object.dataset_names

    @property
    def cut_depth(self):
        return self._cut_depth

    @property
    def hierarchical_clustering_object(self):
        return self._hierarchical_clustering_object

    @property
    def cluster_sizes(self):
        return pd.Series(map(len, self.clusters))

    def __repr__(self):
        if len(self.clusters) <= 10:
            clusters_repr = '\n'.join(map(repr, self.clusters))
        else:
            clusters_repr = '\n'.join(map(repr, self.clusters[:3]))
            clusters_repr += '\n...\n'
            clusters_repr += '\n'.join(map(repr, self.clusters[-3:]))

        return '<ClusterAssignments n={0}, cut depth: {1}\nClusters: \n{2}>'.format(self.n, self.cut_depth, clusters_repr)

    def flatten(self):
        """
        Flattens the data into a `pd.Series` object that gives a cluster number to every element in original data.
        :return:
        """
        buffer = np.empty(self._hierarchical_clustering_object.num_obs)

        for i, cluster in enumerate(self.clusters):
            queue = set([cluster.root])

            while queue:
                node = queue.pop()

                if node.is_leaf():
                    buffer[node.id] = i + 1
                else:
                    queue.add(node.get_left())
                    queue.add(node.get_right())

        return pd.Series(buffer, index=self.hierarchical_clustering_object.data.index)

    @property
    def clusters(self):
        """
        Returns the data clusters in decreasing number of elements
        :return:
        """
        return self._clusters

    def __iter__(self):
        return iter(self.clusters)

    def __getitem__(self, key):
        return self.clusters[key]
