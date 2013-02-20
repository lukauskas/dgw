__author__ = 'saulius'

from multiprocessing import cpu_count, Array, Process, Queue
from Queue import Empty
from math import factorial
import itertools
import numpy as np
import ctypes
from logging import debug

from dgw.dtw.distance import dtw_std

__all__ = ['parallel_pdist']

def combinations_count(n_items):
    """
    Returns the number of distinct combinations of n_items there can be.
    In other words, calculates `len(list(itertools.combinations(n_items,2)))` efficiently.
    :param n_items:
    :return:
    """
    return factorial(n_items) / (2 * factorial(n_items - 2))

def _parallel_pdist_worker(data_buffer, data_buffer_shape, result_buffer, scheduling_queue, exception_queue, dtw_args, dtw_kwargs):
    """
    A worker function for parallel_pdist that is executed on a separate process.

    The function takes a shared memory buffer to read the data from.
    A shared memory buffer to save the data into, `result_buffer` is also passed to the function.

    The function takes the regions of the result function to compute from the scheduling queue `scheduling_queue` and performs
    these calculations until the scheduling function becomes empty.

    Any exceptions that occur are pushed to `exception_queue` and the process exits.

    :param data_buffer: shared memory buffer where data is read from
    :param data_buffer_shape: shape of this buffer -- used to convert it to np array
    :param result_buffer: shared memory buffer to store result in
    :param scheduling_queue: multiprocessing-safe queue to read processing schedule from
    :type scheduling_queue: multiprocessing.Queue
    :param exception_queue: queue that any exceptions that occur will be pushed into
    :param dtw_args: args passed into `dtw_std`
    :param dtw_kwargs: kwargs passed into `dtw_std`
    :return:
    """
    import sys
    import traceback
    import os

    try:
        pid = os.getpid()
        debug('PROCESS {0}: Spawned'.format(pid))

        data_view = np.frombuffer(data_buffer).reshape(data_buffer_shape)  # Point numpy array to memory
        while True:
            schedule = scheduling_queue.get()
            if schedule is None:
                break  # This indicates there is nothing else will be put in the queue. Kind of a "stop" marker
            else:
                debug('PROCESS {0}: Iteration start'.format(pid))
                start, end = schedule

            # Create the combinations object inside the Process so we can just pass start/end locations in the queue
            # Will need to recreate object every time as we use absolute positions to slice
            combs = itertools.combinations(xrange(data_buffer_shape[0]), 2)
            data_indices = itertools.islice(combs, start, end)

            # Actual data processing work
            for i, (x, y) in enumerate(data_indices):
                a = data_view[x]
                b = data_view[y]
                result = dtw_std(a, b, *dtw_args, **dtw_kwargs)
                result_buffer[start + i] = result

            debug('PROCESS {0}: Iteration end'.format(pid))


        debug('PROCESS {0}: Work complete'.format(pid))

    except Exception as e:
        exception_queue.put(e)
        traceback.print_exc()
    except:
        e = sys.exc_info()[0]
        exception_queue.put(e)
        traceback.print_exc()


def parallel_pdist(three_dim_array, n_processes=None, *dtw_args, **dtw_kwargs):
    """
    Calculates pairwise DTW distance for all the rows in three_dim_array provided.
    This module is similar to scipy.spatial.distance.pdist, but uses all CPU cores available, rather than one.
    :param three_dim_array: numpy data array [observations x max(sequence_lengths) x ndim ]
    :param n_processes: number of processes to spawn usage to the number specified.
                        Will default to the number of (virtual) CPUs available if not set
    :param dtw_args: `args` to be passed into `dtw_std`
    :param dtw_kwargs: `kwargs` to be passed into `dtw_std`
    :return: condensed distance matrix (just as `scipy.spatial.distance.pdist`)
    """

    three_dim_array = np.asarray(three_dim_array)
    if n_processes is None:
        n_processes = cpu_count()
    else:
        n_processes = int(n_processes)
        if n_processes <= 0:
            raise ValueError('N_processes should be > 0')
        elif n_processes > cpu_count():
            raise ValueError('The specified number of CPUs to use, {0} is greater than the number of available CPUs, {1}'
                            .format(n_processes, cpu_count()))

    debug('Using {0} processes for parallel computation'.format(n_processes))
    n_items = three_dim_array.shape[0]
    number_of_combinations = combinations_count(n_items)

    # Create a shared memory array to store result
    # Do not lock it as the worker should make sure processes do not overlap the data
    result_buffer = Array(ctypes.c_double, number_of_combinations, lock=False)

    # Create a shared memory buffer for the data array
    shape = three_dim_array.shape
    data_buffer = Array(ctypes.c_double, np.product(shape), lock=False)  # Allocate memory
    data_view = np.frombuffer(data_buffer).reshape(shape)  # Point numpy array to memory
    data_view[:] = three_dim_array  # Copy the contents into the new memory location

    number_of_slices = n_processes * 4

    buffer_size, remainder = divmod(number_of_combinations, number_of_slices)

    scheduling_queue = Queue()

    # Split the data in slices
    for i in xrange(number_of_slices):
        start = i * buffer_size
        end = start + buffer_size
        scheduling_queue.put((start, end))
    if remainder:
        last_start = number_of_slices * buffer_size
        last_end = number_of_combinations
        scheduling_queue.put((last_start, last_end))

    for i in xrange(n_processes):
        scheduling_queue.put(None)  # Add stop items to the queue so we know when its empty for sure

    # Create queue for exceptions
    exception_queue = Queue()

    processes = []
    for i in xrange(n_processes):
        p = Process(target=_parallel_pdist_worker,
                    args=(data_buffer, shape, result_buffer, scheduling_queue, exception_queue, dtw_args, dtw_kwargs))
        processes.append(p)

    # Start all processes
    for p in processes:
        p.start()

    # Join all processes
    for p in processes:
        p.join()

    try:
        exception = exception_queue.get(block=False)
        raise exception
    except Empty:
        pass

    # Convert the result to numpy array in the end
    return np.frombuffer(result_buffer)
