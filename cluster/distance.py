'''
Created on 13 Nov 2012

@author: saulius
'''
from scipy.spatial.distance import pdist
from mlpy import dtw_std

def dtw_distance_matrix(peaks):
    peaks = list(peaks) # Convert to list as need to loop twice
    
    return pdist(peaks, lambda: 0)

