from logging import debug
import fastcluster
import scipy.cluster.hierarchy as hierarchy
from scipy.spatial.distance import num_obs_y
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ..data.containers import AlignmentsData
from ..dtw.distance import dtw_std, no_nans_len
from ..dtw import transformations, dtw_projection

# Inheritting from object explicitly as hierarchy.ClusterNode is not doing this
class DTWClusterNode(object, hierarchy.ClusterNode):
    _prototype = None
    _index = None
    _hierarchical_clustering_object = None
    _projected_data = None
    _warping_paths = None
    _tracked_points = None
    _points_of_interest = None

    def __init__(self, hierarchical_clustering_object, id, prototype, left=None, right=None, dist=0, count=1):
        hierarchy.ClusterNode.__init__(self, id, left=left, right=right, dist=dist, count=count)
        if not isinstance(prototype, pd.DataFrame):
            prototype = pd.DataFrame(prototype, columns=hierarchical_clustering_object.data.dataset_axis)
        self._prototype = prototype
        self._hierarchical_clustering_object = hierarchical_clustering_object
        self._index = pd.Index(self.__get_item_ids())

        # Assume only bin 40 is interesting at this point ... TODO: make this customisable
        self._points_of_interest = {ix: np.array([40], dtype=np.int32) for ix in self.index}

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

    def ensure_points_of_interest_are_tracked_down(self):
        if self._tracked_points is None:
            self._tracked_points = self.__track_points_of_interest()


    @property
    def warping_paths(self):
        if not self._warping_paths:
            self._warping_paths = self.__compute_dtw_warping_paths()
        return self._warping_paths


    def __compute_dtw_warping_paths(self):
        data = self.data
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

    def __project_items_onto_prototype(self):
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

        return AlignmentsData(pd.Panel(projections))

    def __track_points_of_interest(self):
        points_of_interest = self.points_of_interest
        warping_paths = self.warping_paths
        tracked_points = {}
        for ix, poi in points_of_interest.iteritems():
            mapped_points = transformations.points_mapped_to(poi, warping_paths[ix])
            tracked_points[ix] = mapped_points

        return tracked_points
    @property
    def points_of_interest(self):
        poi = self._points_of_interest

        return poi

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

def _compute_dtw_tree(hierarchical_clustering_object, tree, values, prototyping_function):
    """
    Computes prototypes for all nodes in the tree and returns `ClusterNodeWithPrototype` tree
    :param hierarchical_clustering_object: back reference to hierarchical_clustering_object
    :param tree: root of the tree to process
    :param values: values of items in the tree. Should support `values[node.id]`.
    :param prototyping_function: function that will calculate prototypes.
        Must take four arguments: sequence_a, sequence_b, weight_a, weight_b
    :return:
    """
    def map_function(node):
        prototype = values[node.id]

        return DTWClusterNode(hierarchical_clustering_object, node.id, prototype,
                                        left=node.get_left(), right=node.get_right(),
                                        dist=node.dist, count=node.count)

    def reduce_function(parent_node, left, right):
        # Check if we need to swap nodes by chance
        if parent_node.get_left().id == right.id or parent_node.get_right().id == left.id:
            left, right = right, left

        # Consistency check
        assert(parent_node.get_left().id == left.id)
        assert(parent_node.get_right().id == right.id)
        assert(parent_node.count == left.count + right.count)

        prototype = prototyping_function(left.prototype.values, right.prototype.values, left.count, right.count)

        return DTWClusterNode(hierarchical_clustering_object, parent_node.id, prototype,
                                        left=left, right=right, dist=parent_node.dist, count=parent_node.count)

    def is_value(node):
        return isinstance(node, DTWClusterNode)

    return _reduce_tree(tree, reduce_function, map_function, is_value)

class InteractiveDendrogramCutter(object):
    _linkage = None
    _figure = None
    _axis = None

    _cut_value = None
    _line = None

    def _onmotion_listener(self, event):
        if event.inaxes != self._axis or event.ydata is None:
            if self._line.get_visible():
                # Hide line
                self._line.set_visible(False)
                self._figure.canvas.draw()
                return
        else:
            self._line.set_ydata(event.ydata)
            self._line.set_visible(True)
            self._figure.canvas.draw()

    def _onclick_listener(self, event):
        # Allow only double-click, only in axis, only the left button and only regions with value
        if not event.dblclick or event.inaxes != self._axis or event.button != 1 or event.ydata is None:
            return

        self._cut_value = event.ydata
        plt.close(self._figure)

    def __init__(self, linkage):
        self._linkage = linkage


    def show(self):

        figure = plt.figure()
        self._figure = figure
        self._axis = plt.gca()  # Get axis

        self._figure.canvas.mpl_connect('button_press_event', self._onclick_listener)
        self._figure.canvas.mpl_connect('motion_notify_event', self._onmotion_listener)
        hierarchy.dendrogram(self._linkage, no_labels=True)
        self._line = self._axis.axhline(linestyle='--', y=0, visible=False)

        plt.title('Double click on the dendrogram to cut it')
        plt.show()

    @property
    def value(self):
        return self._cut_value

class HierarchicalClustering(object):
    _condensed_distance_matrix = None
    _data = None
    _linkage_matrix = None
    _distance_threshold = None
    __dtw_args = None
    __dtw_kwargs = None

    __tree = None

    def __init__(self, data, condensed_distance_matrix, linkage_matrix=None, dtw_function=dtw_std, prototyping_method='psa'):
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


        .. [#Niennattrakul:2009ep] Vit Niennattrakul and Chotirat Ann Ratanamahatana "Shape averaging under Time Warping",
           2009 6th International Conference on Electrical Engineering/Electronics, Computer,
           Telecommunications and Information Technology (ECTI-CON)

        :param data: a pd.DataFrame object of the data in clusters
        :type data: AlignmentsData
        :param condensed_distance_matrix: a condensed distance matrix of this clustering computed by pdist
        :param linkage_matrix: (optional) linkage_matrix computed by fastcluster.linkage. Will be computed automatically if not provided
        :param dtw_function: DTW calculation function
        :param prototyping_method: Averaging method either 'psa', 'standard' or 'standard-unweighted'
        :return:
        """
        self._condensed_distance_matrix = condensed_distance_matrix
        if not isinstance(data, AlignmentsData):
            raise ValueError('Data should be instance of {0}'.format(AlignmentsData.__name__))

        if not num_obs_y(condensed_distance_matrix) == data.number_of_items:
            raise ValueError('Mismatch of the number of observations '
                             'in the condensed distance matrix and the data provided: '
                             '{0} != {1}'.format(num_obs_y(condensed_distance_matrix), len(data)))

        self._data = data

        # Compute linkage matrix if its not provided
        if linkage_matrix is None:
            linkage_matrix = fastcluster.complete(condensed_distance_matrix)
            # small negative distances in linkage matrix are sometimes possible due to rounding errors. Change them to 0
            linkage_matrix[:, 2][linkage_matrix[:, 2] < 0] = 0
        self._linkage_matrix = linkage_matrix

        self.__dtw_function = dtw_function

        tree = self.__dtw_tree_from_linkage(linkage_matrix, prototyping_method)
        self.__tree = tree

    def __dtw_tree_from_linkage(self, linkage, method):
        """
        Computes a prototyped tree from linkage matrix
        :param linkage: linkage matrix
        :param method: prototyping method
        :return:
        """
        root = hierarchy.to_tree(linkage)
        self._rename_nodes(root)

        if method == 'psa':
            averaging_func = lambda x, y, wx, wy: \
                transformations.sdtw_averaging(x, y, wx, wy, dtw_function=self.dtw_function)
        elif method == 'standard':
            averaging_func = lambda x, y, wx, wy: \
                transformations.dtw_path_averaging(x, y, wx, wy, dtw_function=self.dtw_function)
        elif method == 'standard-unweighted':
            averaging_func = lambda x, y, wx, wy: \
                transformations.dtw_path_averaging(x, y, 1, 1, dtw_function=self.dtw_function)
        else:
            raise ValueError('Incorrect method supplied: '
                             'only \'psa\', \'standard\' or \'standard-unweighted\' supported')

        # Compute prototypes for the tree
        tree = _compute_dtw_tree(self, root, self.data, averaging_func)

        return tree

    @property
    def data(self):
        """
        Returns the data that is internal to the object
        :rtype: AlignmentsData
        """
        return self._data

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

    @property
    def num_obs(self):
        return self.data.number_of_items

    @property
    def dataset_names(self):
        return self.data.dataset_axis

    def as_tree(self):
        return self.__tree

    def dendrogram(self, *args, **kwargs):
        """
        Plots the dendrogram for the hieararchical clustering
        :return:
        """
        linkage = self.linkage
        no_labels = kwargs.pop('no_labels', True)

        color_threshold = kwargs.pop('color_threshold', self._distance_threshold)
        return hierarchy.dendrogram(linkage, no_labels=no_labels, color_threshold=color_threshold, *args, **kwargs)

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
            if current_node.dist > t:
                queue.add(current_node.get_left())
                queue.add(current_node.get_right())
            else:
                clusters.add(current_node)
        return ClusterAssignments(self, clusters, t)

    def interactive_cut(self):
        cutter = InteractiveDendrogramCutter(self.linkage)

        cutter.show()
        value = cutter.value

        try:
            value = float(value)
        except TypeError, ValueError:
            raise ValueError('Incorrect back from the interactive cut routine. Did you double-click it?')

        return self.cut(value)

    def pairwise_distances_to_index(self, query_index):
        """
        Reads the pairwise distances of an index to other indices from the condensed distance matrix
        :param index:
        :return:
        """
        data_index = self.data.items
        if query_index not in data_index:
            raise ValueError('No index {0} in data'.format(query_index))

        n = len(data_index)
        # Find the query index in the full index
        query_index_pos = list(data_index).index(query_index)

        dm = self.condensed_distance_matrix

        dm_indices = np.empty(n-1, dtype=int)
        start = 0
        for i in xrange(query_index_pos):
            # Each i will have to be compared with n_compared to items:
            n_compared_to = n - i - 1

            # the comparisons will be in this order
            # (i, i+1), (i, i+2), (i, i+3)
            # So if index == i+1 then offset = 0
            # If index == i+2, offset = 1
            # So on, so offset = index - i - 1
            query_index_offset = query_index_pos - i-1
            dm_indices[i] = start + query_index_offset

            start += n_compared_to

        # Add all other items
        n_compared_to = n - query_index_pos - 1
        dm_indices[query_index_pos:n-1] = np.arange(start,start+n_compared_to)
        distances = dm[dm_indices]

        data_without_query_index = data_index - [query_index]
        return pd.Series(distances, index=data_without_query_index)



class ClusterAssignments(object):
    _hierarchical_clustering_object = None
    _cluster_roots = None

    _clusters = None
    _cut_depth = None

    def __init__(self, hierarchical_clustering_object, cluster_roots, cut_depth):
        self._hierarchical_clustering_object = hierarchical_clustering_object
        self._cluster_roots = cluster_roots

        clusters = []
        # Store clusters in decreasing number of elements
        for cluster_root in sorted(self._cluster_roots, key=lambda x: x.count, reverse=True):
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
