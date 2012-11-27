'''
Created on 14 Nov 2012

@author: saulius
'''
from scipy.cluster.hierarchy import *
from matplotlib import pyplot

def hierarchical_cluster(data, distance_matrix_func, threshold):
    print "Calculating distance matrix"
    dm = distance_matrix_func(data)
    
    print "Calculating linkage"
    z = linkage(dm, method='complete', metric='cityblock')
    real_threshold = threshold * max(z[:,2])
    dendrogram(z, color_threshold=real_threshold)
    pyplot.savefig('dendrogram.png')
    print "Flattening the clusters"
    cluster_assignments = fcluster(z, t=real_threshold, criterion='distance')
    return cluster_assignments