import fastcluster
import scipy.cluster.hierarchy as hierarchy
from scipy.spatial.distance import num_obs_y
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ..data.containers import AlignmentsData
from ..dtw.distance import dtw_std, no_nans_len
from ..dtw import transformations

def _reduce_tree(tree, reduce_func, map_function=lambda node: node.id):
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
        return reduce_func(left, right)
    ```
    :param tree: root node of the tree
    :type tree: `scipy.cluster.hierarchy.ClusterNode`
    :param reduce_func: function that will be performed on reduction
    :param map_function: function that gets the value of a `ClusterNode` defaults to retrieving the node's id
    :return: The result of reduction function
    """

    def _add_children_to_stack(node):
        stack.append(node.get_left())
        stack.append(node.get_right())

    stack = [tree]

    while stack:
        node = stack.pop()
        if isinstance(node, hierarchy.ClusterNode):
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

            if isinstance(next_node, hierarchy.ClusterNode):
                stack.append(node)
                stack.append(next_node)
            else:
                try:
                    reduced_value = reduce_func(node, next_node)
                except Exception, e:
                    print "Got {0!r} on reduce of {1!r} {2!r}".format(e, node, next_node)
                    raise

                stack.append(reduced_value)

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
        self._axis = plt.gca() # Get axis

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

    def __init__(self, data, condensed_distance_matrix, linkage_matrix=None, dtw_function=dtw_std):
        """
        Initialises hierarchical clustering analyser.
        :param data: a pd.DataFrame object of the data in clusters
        :type data: AlignmentsData
        :param condensed_distance_matrix: a condensed distance matrix of this clustering computed by pdist
        :param linkage_matrix: (optional) linkage_matrix computed by fastcluster.linkage. Will be computed automatically if not provided
        :param dtw_args: arguments that were passed into dtw_std
        :param dtw_kwargs: keyword arguments that were passed into dtw_std
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

        self._linkage_matrix = linkage_matrix

        self.__dtw_function = dtw_function

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

    def as_tree(self):
        root = hierarchy.to_tree(self.linkage)
        self._rename_nodes(root)

        return root

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
            clusters.append(Cluster(hierarchical_clustering_object, cluster_root))
        self._clusters = clusters
        self._cut_depth = cut_depth

    @property
    def n(self):
        return len(self._cluster_roots)
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

class Cluster(object):
    _hierarchical_clustering_object = None
    _index = None
    _item_ids = None

    _cluster_root = None

    _item_dms = None

    def __init__(self, hierarchical_clustering_object, cluster_root):
        """
        Initialises the cluster from the Hiearachical Clustering Object and the root of cluster tree.

        WARNING: The class assumes that the ids of the nodes are the same as indices of the data.
                 See `HierarchicalClustering._rename_nodes` for how this is implemented

        :param hierarchical_clustering_object:
        :type hierarchical_clustering_object: HierarchicalClustering
        :param cluster_root: the root of this cluster
        :type cluster_root: `hierarchy.ClusterNode`
        :return:
        """
        self._hierarchical_clustering_object = hierarchical_clustering_object
        self._cluster_root = cluster_root

        self._item_ids = self.__get_item_ids()
        self._index = pd.Index(self._item_ids)

        # Assert that we did not make up any names that were not in the original index
        assert(len(self._index) == len(hierarchical_clustering_object.data.items & self._index))

    @property
    def hierarchical_clustering_object(self):
        return self._hierarchical_clustering_object

    @property
    def root(self):
        """
        Returns the root of the tree
        :rtype: `hierarchy.ClusterNode`
        """
        return self._cluster_root

    def __get_item_ids(self):
        """
        Gets the ids of all the leaf nodes in the tree.
        :return:
        """

        ids = []
        queue = set([self._cluster_root])

        counter = 0
        while queue:
            node = queue.pop()

            if node.is_leaf():
                ids.append(node.id)
                counter += 1
            else:
                queue.add(node.get_left())
                queue.add(node.get_right())

        return ids

    @property
    def index(self):
        return self._index

    @property
    def n_items(self):
        return len(self.index)

    def __len__(self):
        return self.n_items

    @property
    def items(self):
        """
        Returns a pd.Dataframe of the items in the cluster
        :return:
        """
        return self.hierarchical_clustering_object.data.ix[self.index]

    @property
    def distances(self):
        """
        Returns a DataFrame that lists the distance between each of the clusters
        :return:
        """
        if self._item_dms is None:
            pairwise_distances = map(self.hierarchical_clustering_object.pairwise_distances_to_index, self.index)
            pairwise_distances = map(lambda x: x[self.index], pairwise_distances)
            self._item_dms = pd.DataFrame(pairwise_distances, index=self.index)


        return self._item_dms

    def prototype(self):
        means = self.distances.T.sum()

        min_index_i = np.argmin(means)

        items = self.items
        min_index = items.index[min_index_i]

        return items.ix[min_index]

    def average(self, method='standard'):
        """
        Averages the data in the cluster with respect to hierarchy.

        Supports three methods:

            Prioretised Shape Averaging (PSA) (method='psa')
                as described in [#Niennattrakul:2009ep], see also `dgw.transformations.sdtw_averaging`
            Standard average of DTW path (method='standard')
                averages the DTW paths for each pair of nodes with the same parent.
                The weights are determined by how many sequences were averaged into the node.
                see `dgw.transformations.dtw_path_averaging`
            Unweighted standard average of DTW path (method='standard-unweighted')
                similar to the standard method above, but does not bias sequences higher up the tree


        .. [#Niennattrakul:2009ep] Vit Niennattrakul and Chotirat Ann Ratanamahatana "Shape averaging under Time Warping",
           2009 6th International Conference on Electrical Engineering/Electronics, Computer,
           Telecommunications and Information Technology (ECTI-CON)

        :param method: method to use either 'psa' or 'standard'
        :return:
        """
        if method == 'psa':
            average_item = self.__average_sdtw()
        elif method == 'standard':
            average_item = self.__average_standard()
        elif method == 'standard-unweighted':
            average_item = self._average_standard_unweighted()
        else:
            raise ValueError('Incorrect method supplied: only \'psa\', \'standard\' or \'standard-unweighted\' supported')

        df = pd.DataFrame(average_item, columns=self.hierarchical_clustering_object.data.dataset_axis)
        return df

    def _average_standard_unweighted(self):

        def reduce_function(x, y):
            path_average = transformations.dtw_path_averaging(x, y, dtw_function=self.dtw_function)
            path_average = transformations.uniform_shrinking_to_length(path_average,
                                                                               max(no_nans_len(x), no_nans_len(y)))
            return path_average

        items_ix = self.items.ix
        map_function = lambda x: items_ix[x.id].values

        final_item = _reduce_tree(self.root, reduce_function, map_function)

        return final_item
    @property
    def dtw_function(self):
        return self.hierarchical_clustering_object.dtw_function

    def __average_standard(self):
        def reduce_function(x, y):
            sequence_a, weight_a = x
            sequence_b, weight_b = y
            path_average = transformations.dtw_path_averaging(sequence_a, sequence_b, weight_a, weight_b,
                                                                      dtw_function=self.dtw_function)

            path_average = transformations.uniform_shrinking_to_length(path_average,
                                                                               max(no_nans_len(sequence_a),
                                                                                   no_nans_len(sequence_b)))
            return path_average, \
                   weight_a + weight_b

        return self.__weighted_averaging(reduce_function)

    def __average_sdtw(self):
        # This is essentially Prioritised Shape Averaging as described in Vit Niennattrakul and Chotirat Ann Ratanamahatana
        def reduce_function(x, y):
            sequence_a, weight_a = x
            sequence_b, weight_b = y
            avg_seq = transformations.sdtw_averaging(sequence_a, sequence_b, weight_a, weight_b,
                                                             dtw_function=self.dtw_function)
            return avg_seq, weight_a + weight_b

        return self.__weighted_averaging(reduce_function)

    def __weighted_averaging(self, reduce_function):
        items_ix = self.items.ix
        map_function = lambda x: (items_ix[x.id].values, 1)

        final_item = _reduce_tree(self.root, reduce_function, map_function)
        ans = final_item[0]

        return ans

    def __repr__(self):
        return '<Cluster n_items={0}>'.format(self.n_items)
