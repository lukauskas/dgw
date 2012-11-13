from __future__ import print_function
import matplotlib.pyplot as plt

import numpy as np

def draw_heatmap(adjusted_peak_data):
    adjusted_peak_data = list(adjusted_peak_data)
    adjusted_peak_data = map(dict, adjusted_peak_data)
    
    min_offset = None
    max_offset = None
    
    min_value = None
    max_value = None
    
    for apd in adjusted_peak_data:
        for offset, val in apd.iteritems():
            if offset < min_offset or min_offset is None:
                min_offset = offset
            
            if offset > max_offset or max_offset is None:
                max_offset = offset
                
            if val < min_value or min_value is None:
                min_value = val
            
            if val > max_value or max_value is None:
                max_value = val
            
            
            
    N = len(adjusted_peak_data)

    arr = np.empty([N, max_offset+1-min_offset], dtype=float)
    
    for i, apd in enumerate(adjusted_peak_data):
        apd = dict(apd)
        
        for (pos, n) in sorted(apd.iteritems()):
            pos = pos - min_offset
            scaled_n = (float(n) - min_value) / (max_value - min_value)
            
            arr[i][pos] = scaled_n
    plt.imshow(arr)
    plt.show()