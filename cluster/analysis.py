__author__ = 'saulius'
import fastcluster
import scipy.cluster.hierarchy as hierarchy
from scipy.spatial.distance import num_obs_y
import pandas as pd


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
        cluster_assignments = hierarchy.fcluster(self.linkage, t=t, criterion=criterion, *args, **kwargs)
        clusters = pd.Series(cluster_assignments, index=self._data.index)
        return ClusterAssignments(self, clusters)

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
        ans = []
        cluster_assignments = self._cluster_assignments
        hco = self.hierarchical_clustering_object
        for cluster_i in cluster_assignments:
            indices = cluster_assignments[cluster_assignments==cluster_i].index
            ans.append(Cluster(hco, indices))

        return ans

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
        return self.hierarchical_clustering_object.data.ix[self.item_indices]

    def __repr__(self):
        '<Cluster n_items={0} >'.format(self.n_items)
