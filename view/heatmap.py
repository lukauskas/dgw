from __future__ import print_function
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib
import pandas as pd

def plot(data_frame):

    # Calculate the max_value threshold so heatmap looks good
    values = data_frame.unstack().dropna()
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

    # Start plotting
    N = len(data_frame)

    # Create a masked array from data_frame
    masked_arr = np.ma.array(data_frame, mask=np.isnan(data_frame))

    # Plot the array using JET colourmap, draw NaNs as grey
    cmap = matplotlib.cm.get_cmap('jet')
    cmap.set_bad('#CCCCCC', 1.0)
    
    result = plt.pcolormesh(masked_arr, cmap=cmap)
    result.set_clim(values.min(), max_value)
    plt.colorbar(result)
    ax = plt.gca()

#    labels = data_frame.columns
#    ax.set_xticklabels(labels)

    #plt.xticks(plt.xticks, map(lambda x: x-min_offset, plt.xticks))
    #plt.savefig("heatmap.png", bbox_inches='tight')
    #plt.show()