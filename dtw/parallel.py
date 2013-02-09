__author__ = 'saulius'

from distance import dtw_std
from multiprocessing import Pool
from math import factorial
import itertools
import numpy as np
import gc

class combinations_with_len(itertools.combinations):
    """
    Extends itertools.combinations by providing __len__ attribute.
    This allows better integration with multiprocessing.Pool
    but it means that all iterables passed to __init__ method must have the __len__ attribute as well
    """
    __length = None

    def __calculate_length(self, iterable_len):
        return factorial(iterable_len) / (2 * factorial(iterable_len - 2))

    def __init__(self, iterable_with_len, r):
        super(combinations_with_len, self).__init__(iterable_with_len, r)
        iterable_length = len(iterable_with_len) # Should raise TypeError if iterable does not have __len__
        self.__length = self.__calculate_length(iterable_length)

    def __len__(self):
        return self.__length

def parallel_pdist(three_dim_array, metric='sqeuclidean'):
    """
    Calculates pairwise DTW distance for all the rows in three_dim_array provided.
    This module is similar to scipy.spatial.distance.pdist, but uses all CPU cores available, rather than one.
    :param three_dim_array: numpy data array [observations x max(sequence_lengths) x ndim ]
    :param metric: either 'euclidean' or 'sqeuclidean'
    :return: condensed distance matrix (just as pdist)
    """

    three_dim_array = np.asarray(three_dim_array)
    p = Pool()
    smd = _shared_mem_dtw(three_dim_array, metric=metric)

    n_items = three_dim_array.shape[0]
    combs = combinations_with_len(range(n_items), 2)
    ans = p.map(smd, combs)

    p.close()
    p.join()
    return np.array(ans)

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