'''
Created on 14 Nov 2012

@author: saulius
'''
from scipy.cluster.hierarchy import *
from matplotlib import pyplot

def hierarchical_cluster(data, distance_matrix_func, threshold, z = None, **kwargs):
    print "Calculating distance matrix"
    dm = distance_matrix_func(data)
    
    print "Calculating linkage"
    if z is None:
        z = linkage(dm, method='complete', metric='euclidean')
    
    y = inconsistent(z, 3)
    f = open('inconsistency.txt', 'w')
    for something in y:
        f.write('{0}\n'.format('\t'.join(map(str, something))))
    f.close()
    dendrogram(z, color_threshold=threshold)
    pyplot.savefig('dendrogram.png')
    print "Flattening the clusters"
    cluster_assignments = fcluster(z, t=threshold, **kwargs)
    return cluster_assignments, z