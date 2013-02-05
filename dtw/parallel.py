__author__ = 'saulius'

from distance import dtw_std
from multiprocessing import Pool
from math import factorial
from itertools import combinations
import numpy as np
import gc

def parallel_pdist(two_dim_array, metric='sqeuclidean'):
    '''
    Calculates pairwise DTW distance for all the rows in
    two_dim_array using all CPUs of the computer.

    Note: if using pandas dataframes, call the function as follows:
        parallel_pdist(df.values)
    otherwise some weird behaviour will happen.

    @param two_dim_array: a numpy data array.
    @return: condensed distance matrix. See pdist documentation in scipy
    '''
    two_dim_array = np.asarray(two_dim_array)
    p = Pool()
    smd = _shared_mem_dtw(two_dim_array, metric=metric)

    size_dm = factorial(len(two_dim_array)) / (2 * factorial(len(two_dim_array) - 2))
    combs = combinations(xrange(len(two_dim_array)), 2)

    step_size = 30000000 # TODO: this should be related to total/free memory somehow

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

class _shared_mem_dtw:
    '''
        Dtw function that keeps the data matrix in memory
        and takes only indexes as call arguments
    '''
    shared_mem_matrix = None
    metric = None
    def __init__(self, shared_mem_matrix, metric='sqeuclidean'):

        self.metric = metric
        self.shared_mem_matrix = shared_mem_matrix

    def __call__(self, args):

        x = args[0]
        y = args[1]

        smm = self.shared_mem_matrix

        a = smm[x]
        b = smm[y]

        return dtw_std(a,b, metric=self.metric)