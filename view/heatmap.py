from __future__ import print_function
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib

def draw_heatmap(adjusted_peak_data, ax=None):
    adjusted_peak_data = list(adjusted_peak_data)
    min_offset = None
    max_offset = None
    
    min_value = None
    max_value = None

    values = []
    for apd in adjusted_peak_data:
        for offset, val in apd:
            if offset < min_offset or min_offset is None:
                min_offset = offset
            
            if offset > max_offset or max_offset is None:
                max_offset = offset
                
            if val < min_value or min_value is None:
                min_value = val
            
            if val > max_value or max_value is None:
                max_value = val
            
            values.append(val)
   
    min_value = min(values)
    max_value = max(values)
    
    print("calculating histogram")
    histogram, bin_edges = np.histogram(values, bins=1000)
    last_good_bin = None
    TRIM_MAX_VALUE_THRESHOLD = 0.001 # Trim max value threshold to be the value that hides 0.1% of peaks
    allowed_slack = round(TRIM_MAX_VALUE_THRESHOLD * len(values))
    current_slack = 0
    for i, hist in reversed(list(enumerate(histogram))):
        current_slack =+ hist
        if current_slack >= allowed_slack:
            last_good_bin = i
            break
    
    max_value = bin_edges[last_good_bin+1]
          
    N = len(adjusted_peak_data)

    # Initialise empty array that will later become our heatmap
    arr = np.empty([N, max_offset+1-min_offset], dtype=float)
    arr[:] = np.NAN
    
    # Fill in the values
    # Sort adjusted peak data by the length of peaks
    adjusted_peak_data = sorted(adjusted_peak_data, key = len)
    for i, apd in enumerate(adjusted_peak_data):
        
        for (pos, n) in apd:
            pos = pos - min_offset
            #scaled_n = (float(n) - min_value) / (max_value - min_value)
            
            arr[i,pos] = n
    

    # Create a masked array from it
    masked_arr = np.ma.array(arr, mask=np.isnan(arr))
    
    
    # Plot the array using JET Cmap, draw NaNs as black
    cmap = matplotlib.cm.get_cmap('jet')
    cmap.set_bad('#CCCCCC', 1.0)
    
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    
    result = ax.imshow(masked_arr, cmap=cmap, aspect="auto", interpolation="nearest")
    print('Clipping result\'s colours to: ({0}, {1})'.format(min_value, max_value))
    result.set_clim(min_value, max_value)
    cbar = plt.colorbar(result)
    
    #plt.savefig("heatmap.png", bbox_inches='tight')
    #plt.show()