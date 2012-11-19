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


def dtw(x, y):
    pass

def shrink_time_series(x, shrink_factor):
    '''
    Shrinks time series x by a factor specified as shrink_factor.
    The series is shrunk by averaging components.
    
    Examples:
    >>> shrink_time_series([1,2,3,4,5,6,7,8,9], 2.0)
    [1.5, 4.0, 6.5, 8.5]
    >>> shrink_time_series([1,2,3,4,5,6], 3.0)
    [2.0, 5.0]

    This is direct translation from a similar function in PAA.java in FastDTW implementation.
    
    @param x: time series to be shrunk. should be list like and support direct element access
    @param shrink_factor: ratio to reduce the size by. E.g
    @type shrink_factor: float
    '''
    
    assert(shrink_factor > 1)
    
    original_size = len(x)
    # Get the size of shrunk time series
    reduced_size = int(float(original_size) / shrink_factor)
    
    # Determine size of a sampled font (might be fractional)
    reduced_point_size = float(original_size) / reduced_size
   
    # Keep track of points being averaged
    read_from = 0
    read_to   = None
    
    ans = []
    while (read_from < original_size):
        
        read_to = int(round(reduced_point_size*(len(ans)+1)) -1 )
    
        points_to_read = read_to - read_from +1
        
        val_sum = 0
        
        for i in range(read_from, read_to+1):
            val_sum += x[i]
        
        val_sum = float(val_sum) / points_to_read
        ans.append(val_sum)
        read_from = read_to + 1
        
    return ans

def expanded_res_window(low_res_path, x, y, radius):
    pass


def fast_dtw(x, y, radius):
    
    min_ts_size = radius + 2
    
    len_x = len(x)
    len_y = len(y)
    if (len_x <= min_ts_size or len_y <= min_ts_size):
        return dtw(x, y)
    else:
        SHRINK_FACTOR = 2.0
        shrunk_x = shrink_time_series(x, SHRINK_FACTOR)
        shrunk_y = shrink_time_series(x, SHRINK_FACTOR)
        
        low_res_path = fast_dtw(shrunk_x, shrunk_y, radius)
        
        window = expanded_res_window(low_res_path, x, y, radius)

        return dtw(x,y, window)
