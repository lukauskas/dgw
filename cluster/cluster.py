'''
Created on 14 Nov 2012

@author: saulius
'''
from scipy.cluster.hierarchy import *
from matplotlib import pyplot

def hierarchical_cluster(data, distance_matrix_func, threshold):

    dm = distance_matrix_func(data)
    
    z = linkage(dm)
    
    cluster_assignments = fcluster(z, t=threshold, criterion='distance')
    return cluster_assignments