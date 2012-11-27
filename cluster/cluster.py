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
    dendrogram(z)
    pyplot.savefig('dendrogram.png')
    print "Flattening the clusters"
    cluster_assignments = fcluster(z, t=threshold, criterion='maxclust')
    return cluster_assignments