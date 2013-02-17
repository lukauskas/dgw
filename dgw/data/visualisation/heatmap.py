from __future__ import print_function
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import matplotlib
import pandas as pd

def plot(data_frame, sort_by='length', clip_colors=True):
    values = data_frame.unstack().dropna()
    if clip_colors:
        # Calculate the max_value threshold so heatmap looks good

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
    else:
        max_value = values.max()

    lengths = data_frame.apply(lambda x : len(x.values[np.invert(np.isnan(x.values))]), axis=1)
    lengths.sort()

    # Cut most of the columns out of dataframe
    data_frame = data_frame[data_frame.columns[:lengths.max()+1]]
    if sort_by == 'length':
        data_frame = data_frame.ix[lengths.index]
    elif sort_by == 'read_count':
        read_counts = data_frame.T.sum()
        read_counts.sort()

        data_frame = data_frame.ix[read_counts.index]

    elif sort_by == 'nothing':
        pass
    else:
        raise ValueError("Only 'length', 'read_count' and 'nothing' supported")

    # Start plotting
    N = len(data_frame)

    # Create a masked array from data_frame
    masked_arr = np.ma.array(data_frame, mask=np.isnan(data_frame))

    # Plot the array using JET colourmap, draw NaNs as grey
    cmap = matplotlib.cm.get_cmap('jet')
    cmap.set_bad('#FFFFFF', 1.0)
    
    result = plt.imshow(masked_arr, cmap=cmap, aspect="auto", interpolation='nearest')
    result.set_clim(values.min(), max_value)
    plt.colorbar(result)
    ax = plt.gca()

#    labels = data_frame.columns
#    ax.set_xticklabels(labels)

    #plt.xticks(plt.xticks, map(lambda x: x-min_offset, plt.xticks))
    #plt.savefig("heatmap.png", bbox_inches='tight')
    #plt.show()