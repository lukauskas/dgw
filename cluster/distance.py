'''
Created on 13 Nov 2012

@author: saulius
'''
from scipy.spatial.distance import pdist
from mlpy import dtw_std

def dtw_distance_matrix(peaks):
    peaks = list(peaks) # Convert to list as need to loop twice
    
    distances = []
    while (peaks):
        
        head = peaks[0]
        tail = peaks[1:]
        peaks = tail
        
        for peak in tail:
            distances.append(dtw_std(head, peak))
    
    return distances
