import numpy as np
from mlpy import dtw_std as mlpy_dtw_std


import gc
from itertools import combinations
from multiprocessing import Pool as Pool
import numpy as np

from math import factorial

def dtw_std(x, y, *args, **kwargs):
    '''
        Wrapper around mlpy's dtw_std that first strips all NaNs out of the data.

    @param x:
    @param y:
    @param args:
    @param kwargs:
    @return:
    '''
    '''
    @param x:
    @param y:
    @param args:
    @param kwargs:
    @return:
    '''
    x = x[np.invert(np.isnan(x))]
    y = y[np.invert(np.isnan(y))]

    return mlpy_dtw_std(x, y, *args, **kwargs)

#-----------------------------------------------------------------------------------------------------------------------

def parallel_pdist(two_dim_array):
    '''
    Calculates pairwise DTW distance for all the rows in
    two_dim_array using all CPUs of the computer.

    Note: if using pandas dataframes, call the function as follows:
        parallel_pdist(df.values)
    otherwise some weird behaviour will happen.
    
    @param two_dim_array: a numpy data array.
    @return: condensed distance matrix. See pdist documentation in scipy
    '''

    p = Pool()
    smd = _shared_mem_dtw(two_dim_array)

    size_dm = factorial(len(two_dim_array)) / (2 * factorial(len(two_dim_array) - 2))
    combs = combinations(xrange(len(two_dim_array)), 2)

    step_size = 30000000
    ans_container = np.empty(size_dm)
    curr_step_offset = 0
    while True:
        items_to_process = take(combs, step_size)

        items_to_process = list(items_to_process)
        if not items_to_process:
            break

        ans = p.map(smd, items_to_process)

        start = curr_step_offset*step_size
        end   = start + len(items_to_process)

        ans_container[start:end] = ans
        curr_step_offset += 1

        del items_to_process
        gc.collect()

    p.close()
    p.join()
    return ans_container



class _shared_mem_dtw():
    '''
        Dtw function that keeps the data matrix in memory
        and takes only indexes as call arguments
    '''
    shared_mem_matrix = None
    def __init__(self, shared_mem_matrix):
        self.shared_mem_matrix = shared_mem_matrix

    def __call__(self, args):

        x = args[0]
        y = args[1]

        a = self.shared_mem_matrix[x]
        b = self.shared_mem_matrix[y]

        return dtw_std(a,b)

def take(iterable, n):
    '''
    Take first n elements from iterable
    @param iterable:
    @param n:
    @return:
    '''
    counter = 0

    for x in iterable:
        counter += 1
        yield x
        if counter >= n:
            break
