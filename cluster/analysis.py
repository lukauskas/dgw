__author__ = 'saulius'
import fastcluster
import scipy.cluster.hierarchy as hierarchy
from scipy.spatial.distance import num_obs_y
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

    def __init__(self, data, condensed_distance_matrix, linkage_matrix=None):
        """
        Initialises hierarchical clustering analyser.
        :param data: a pd.DataFrame object of the data in clusters
        :param condensed_distance_matrix: a condensed distance matrix of this clustering computed by pdist
        :param linkage_matrix: (optional) linkage_matrix computed by fastcluster.linkage. Will be computed automatically if not provided
        :return:
        """
        self._condensed_distance_matrix = condensed_distance_matrix
        if not isinstance(data, pd.DataFrame):
            raise ValueError('Data should be instance of pd.DataFrame')

        if not num_obs_y(condensed_distance_matrix) == len(data):
            raise ValueError('Mismatch of the number of observations '
                             'in the condensed distance matrix and the data provided: '
                             '{0} != {1}'.format(num_obs_y(condensed_distance_matrix), len(data)))

        self._data = data
        self._linkage_matrix = linkage_matrix

    @property
    def data(self):
        return self._data

    @property
    def condensed_distance_matrix(self):
        return self._condensed_distance_matrix

    @property
    def linkage(self):
        """
        Returns the linkage matrix of the clustering provided
        :return: linkage matrix.
        """
        if self._linkage_matrix is None:
            # If we haven't already, compute complete linkage matrix
            self._linkage_matrix = fastcluster.complete(self._condensed_distance_matrix)

        return self._linkage_matrix

    def dendrogram(self, *args, **kwargs):
        """
        Plots the dendrogram for the hieararchical clustering
        :return:
        """
        linkage = self.linkage
        no_labels = kwargs.pop('no_labels', True)

        return hierarchy.dendrogram(linkage, no_labels=no_labels, *args, **kwargs)

    def cut(self, t, criterion='distance', *args, **kwargs):
        """
        Cuts the dendrogram at specified threshold t.
        See scipy.cluster.hierarchy.fcluster
        :param t: threshold
        :param criterion:
        :param args:
        :param kwargs:
        :return:
        """
        cluster_assignments = hierarchy.fcluster(self.linkage, t=t, criterion=criterion, *args, **kwargs)
        clusters = pd.Series(cluster_assignments, index=self._data.index)
        return ClusterAssignments(self, clusters)

    def interactive_cut(self):
        cutter = InteractiveDendrogramCutter(self.linkage)

        cutter.show()
        value = cutter.value

        try:
            value = float(value)
        except TypeError, ValueError:
            raise ValueError('Incorrect back from the interactive cut routine. Did you double-click it?')

        return self.cut(value, criterion='distance')

    def pairwise_distances_to_index(self, query_index):
        """
        Reads the pairwise distances of an index to other indices from the condensed distance matrix
        :param index:
        :return:
        """
        data_index = self.data.index
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
    _cluster_assignments = None

    _clusters = None

    def __init__(self, hierarchical_clustering_object, cluster_assignments):
        self._hierarchical_clustering_object = hierarchical_clustering_object

        if not isinstance(cluster_assignments, pd.Series):
            cluster_assignments = pd.Series(cluster_assignments, self._hierarchical_clustering_object.data.index)

        self._cluster_assignments = cluster_assignments

    @property
    def n(self):
        return len(self._cluster_assignments.value_counts())

    @property
    def hierarchical_clustering_object(self):
        return self._hierarchical_clustering_object

    def __repr__(self):
        return '<ClusterAssignments n={0} for {1!r}>'.format(self.n, self._hierarchical_clustering_object)

    @property
    def clusters(self):
        """
        Returns the data clusters in decreasing order of number of elements
        :return:
        """

        if self._clusters is None:
            clusters = []
            cluster_assignments = self._cluster_assignments
            hco = self.hierarchical_clustering_object
            for cluster_i in cluster_assignments.value_counts().index:
                indices = cluster_assignments[cluster_assignments==cluster_i].index
                clusters.append(Cluster(hco, indices))
            self._clusters = clusters

        return self._clusters

    def __iter__(self):
        return iter(self.clusters)

    def __getitem__(self, key):
        return self.clusters[key]

class Cluster(object):
    _hierarchical_clustering_object = None
    _item_indices = None

    _item_dms = None

    def __init__(self, hierarchical_clustering_object, item_indices):
        """

        :param hierarchical_clustering_object:
        :type hierarchical_clustering_object: HierarchicalClustering
        :param item_indices:
        :return:
        """
        self._hierarchical_clustering_object = hierarchical_clustering_object
        self._item_indices = item_indices

    @property
    def hierarchical_clustering_object(self):
        return self._hierarchical_clustering_object

    @property
    def index(self):
        return self._item_indices

    @property
    def n_items(self):
        return len(self._item_indices)

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
            pairwise_distances = map( lambda x: x[self.index], pairwise_distances)
            self._item_dms = pd.DataFrame(pairwise_distances, index=self.index)


        return self._item_dms

    def prototype_item(self): # TODO: Consider renaming to prototype
        means = self.distances.T.mean()

        min_index_i = np.argmin(means)

        items = self.items
        min_index = items.index[min_index_i]

        return items.ix[min_index]

    def __repr__(self):
        return '<Cluster n_items={0} >'.format(self.n_items)
