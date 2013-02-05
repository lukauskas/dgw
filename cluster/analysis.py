__author__ = 'saulius'
import fastcluster
import scipy.cluster.hierarchy as hierarchy
from scipy.spatial.distance import num_obs_y
import pandas as pd

import matplotlib.pyplot as plt

class InteractiveDendrogramCutter(object):
    _linkage = None
    _figure = None

    _cut_value = None

    def _onclick_listener(self, event):
        # Allow only double-click, only in axis, only the left button and only regions with value
        if not event.dblclick or event.inaxes is None or event.button != 1 or event.ydata is None:
            return

        self._cut_value = event.ydata
        plt.close(self._figure)

    def __init__(self, linkage, figure=None):
        self._linkage = linkage
        if figure is None:
            figure = plt.figure()
        self._figure = figure

        self._figure.canvas.mpl_connect('button_press_event', self._onclick_listener)

    def show(self):
        hierarchy.dendrogram(self._linkage, no_labels=True)
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

class ClusterAssignments(object):
    _hierarchical_clustering_object = None
    _cluster_assignments = None

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
        ans = []
        cluster_assignments = self._cluster_assignments
        hco = self.hierarchical_clustering_object
        for cluster_i in cluster_assignments.value_counts().index:
            indices = cluster_assignments[cluster_assignments==cluster_i].index
            ans.append(Cluster(hco, indices))

        return ans

    def __iter__(self):
        return iter(self.clusters)

class Cluster(object):
    _hierarchical_clustering_object = None
    _item_indices = None

    def __init__(self, hierarchical_clustering_object, item_indices):
        self._hierarchical_clustering_object = hierarchical_clustering_object
        self._item_indices = item_indices

    @property
    def hierarchical_clustering_object(self):
        return self._hierarchical_clustering_object

    @property
    def item_indices(self):
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
        return self.hierarchical_clustering_object.data.ix[self.item_indices]

    def __repr__(self):
        return '<Cluster n_items={0} >'.format(self.n_items)
