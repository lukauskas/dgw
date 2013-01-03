from collections import defaultdict
import numpy as np
from mlpy import dtw_std as mlpy_dtw_std


import gc
from itertools import combinations, izip
from multiprocessing import Pool as Pool
import numpy as np
import matplotlib.pyplot as plt

from math import factorial

def dtw_std(x, y, metric='sqeuclidean', *args, **kwargs):
    '''
        Wrapper around mlpy's dtw_std that first strips all NaNs out of the data.

    @param x:
    @param y:
    @param args:
    @param kwargs:
    @return:
    '''

    x = _strip_nans(x)
    y = _strip_nans(y)

    if metric == 'sqeuclidean':
        squared = True
    elif metric == 'euclidean':
        squared = False
    else:
        raise ValueError('Unsupported metric provided: {0!r}'.format(metric))

    return mlpy_dtw_std(x, y, squared=squared, *args, **kwargs)

def adaptive_scaling(x, threshold=1e-6):
    '''
        Merges successive coordinates that are identical (i.e. difference is smaller than threshold provided)

        See:
        Petitjean, F., Ketterlin, A., & Gancarski, P. (2011).
        A global averaging method for dynamic time warping, with applications to clustering.
        Pattern Recognition, 44(3), 678-693. doi:10.1016/j.patcog.2010.09.013

    :param x:
    :param threshold:
    :return:
    '''

    new_sequence = []
    prev_item = None
    for item in x:
        if np.isnan(item):
            break
        if prev_item is None or abs(item - prev_item) > threshold:
            new_sequence.append(item)
        prev_item = item

    # Append NaNs to the end.
    ans = np.empty(len(x))
    ans[:len(new_sequence)] = new_sequence
    ans[len(new_sequence):] = np.nan

    return ans

def _strip_nans(sequence):
    '''
        Strips NaNs that are padded to the right of each peak if they are of unequal length
    :param sequence:
    :return:
    '''
    return sequence[np.invert(np.isnan(sequence))]

#-----------------------------------------------------------------------------------------------------------------------
def dba(sequences, initialisation_sequence, convergence_threshold=1e-6, metric='sqeuclidean'):
    '''
        Implements DBA (DTW Barycenter Averaging) algorithm given by
        Petitjean, F., Ketterlin, A., & Gancarski, P. (2011).
        A global averaging method for dynamic time warping, with applications to clustering.
        Pattern Recognition, 44(3), 678-693. doi:10.1016/j.patcog.2010.09.013

    :param sequences: sequences in the cluster
    :param initialisation_sequence: initialisation sequence to choose
    :param convergence_threshold: threshold at which values will be considered equal and iteration will stop
    :param metric: metric to use either euclidean or sqeuclidean
    :return:
    '''

    MAX_ITERATIONS = 2000 # How many iterations of the algorithm to allow

    # Init iteration
    previous_base = initialisation_sequence
    iteration = 1
    new_base = None

    completed = False

    while iteration <= MAX_ITERATIONS:

        total_distance = 0
        base_assoc = {}
        for sequence in sequences:

            distance, cost, path = dtw_std(previous_base, sequence, dist_only=False)
            total_distance += distance

            # splitting and zipping again is necessary here
            # as numpy arrays are not too friendly for iteration like this
            path_base, path_sequence = path
            for (base_coord, sequence_coord) in zip(path_base, path_sequence):
                # We essentially compute barycenter here which is defined as
                # barycenter(x_1, ..., x_n) = SUM(x_1,..,x_n) / n
                # Where x_1 and x_n are given by assoc(S1,..S2) that links
                # that links each coordinate of the average sequence to one or more coordinates of the sequences
                # of S1..S2 see Petitjean et al, 2011

                try:
                    base_assoc[base_coord].append(sequence[sequence_coord])
                except KeyError:
                    base_assoc[base_coord] = [sequence[sequence_coord]]

        #new_base = sums_of_values / counts_of_mappings # For euclidean metric
        # We need a median to minimise unsquared city-block distance
        if metric == 'euclidean':
            new_base = np.array([ np.median(base_assoc[i]) if not np.isnan(previous_base[i]) else np.nan for i in range(len(previous_base))])
        elif metric == 'sqeuclidean':
            new_base = np.array([ np.mean(base_assoc[i]) if not np.isnan(previous_base[i]) else np.nan for i in range(len(previous_base))])
        else:
            raise ValueError('Metric {0!r} is unsupported'.format(metric))

        # If we converged already, return
        difference = np.abs(_strip_nans(new_base) - _strip_nans(previous_base))
        #print iteration, np.sum(difference), np.min(difference), np.max(difference), np.mean(difference)
        if (difference < convergence_threshold).all():
            completed = True
            break
        else:
            # Else, go further down the rabbit hole and continue
            iteration += 1
            previous_base = new_base

    if not completed:
        raise Exception('Maximum limit of {0} iterations has been reached'.format(MAX_ITERATIONS))
    return new_base, total_distance

def dtw_projection(base, x):
    '''
    Projects time series x onto a base time series using dtw
    :param base: base time series
    :param x: another time series to project on
    :return: new time series of length base containing x projected on it
    '''

    distance, cost, path = dtw_std(base, x, dist_only=False)

    path_base, path_other = path

    current_sums = np.zeros(len(base))
    current_counts = np.zeros(len(base))

    current_sums[np.isnan(base)] = np.nan
    current_counts[np.isnan(base)] = np.nan

    for mapped_i, i in zip(path_base, path_other):
        # Go through the path and sum all points that map to the base location i together
        current_sums[mapped_i] += x[i]
        current_counts[mapped_i] += 1


    current_average = current_sums / current_counts

    return current_average

def min_dist_to_others(dm, n):

    sums = defaultdict(lambda: 0)
    for (a, b), dist in izip(combinations(xrange(n), 2), dm):
        sums[a] += dist
        sums[b] += dist

    min_val = float('+inf')
    min_i = None
    for i, val in sums.iteritems():
        if val < min_val:
            min_val = val
            min_i = i

    return min_i, min_val

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

    p = Pool()
    smd = _shared_mem_dtw(two_dim_array, metric=metric)

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

        a = self.shared_mem_matrix[x]
        b = self.shared_mem_matrix[y]

        return dtw_std(a,b, metric=self.metric)

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
