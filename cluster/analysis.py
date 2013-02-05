__author__ = 'saulius'
import fastcluster
import scipy.cluster.hierarchy as hierarchy
from scipy.spatial.distance import num_obs_y
import pandas as pd


class HierarchicalClustering(object):
    dm = None
    data = None
    _linkage_matrix = None

    def __init__(self, data, condensed_distance_matrix, linkage_matrix=None):
        """
        Initialises hierarchical clustering analyser.
        :param data: a pd.DataFrame object of the data in clusters
        :param condensed_distance_matrix: a condensed distance matrix of this clustering computed by pdist
        :param linkage_matrix: (optional) linkage_matrix computed by fastcluster.linkage. Will be computed automatically if not provided
        :return:
        """
        self.dm = condensed_distance_matrix
        if not isinstance(data, pd.DataFrame):
            raise ValueError('Data should be instance of pd.DataFrame')

        if not num_obs_y(condensed_distance_matrix) == len(data):
            raise ValueError('Mismatch of the number of observations '
                             'in the condensed distance matrix and the data provided: '
                             '{0} != {1}'.format(num_obs_y(condensed_distance_matrix), len(data)))

        self.data = data
        self._linkage_matrix = linkage_matrix

    @property
    def linkage(self):
        """
        Returns the linkage matrix of the clustering provided
        :return: linkage matrix.
        """
        if self._linkage_matrix is None:
            # If we haven't already, compute complete linkage matrix
            self._linkage_matrix = fastcluster.complete(self.dm)

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
        clusters = pd.Series(cluster_assignments, index=self.data.index)
        return clusters




