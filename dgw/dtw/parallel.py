from multiprocessing import cpu_count, Array, Process, Queue
from Queue import Empty, Full
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

def _parallel_dtw_worker(data_buffer, operations_generator, data_buffer_shape, result_buffer,
                         scheduling_queue, exception_queue,
                         dtw_args, dtw_kwargs):
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
    :param dtw_args: args passed into `dtw_std` or `uniform_scaled_distance`
    :param dtw_kwargs: kwargs passed into `dtw_std` or `uniform_scaled_distance`
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
            combs = operations_generator(xrange(data_buffer_shape[0]))
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

def _pdist_operations_generator_factory(indices):
    return itertools.combinations(indices, 2)

def _parallel_dtw(three_dim_array, _operations_generator_factory, n_operations, n_processes=None, *dtw_args, **dtw_kwargs):
    """
    Runs DTW on parallel
    :param three_dim_array: three-dimensional numpy array of data
    :param _operations_generator_factory: factory function that generates the operations required.
       Should take one argument -- all indices in the data and return a generator of [(i1, j1), (i1,j2), ...]
       where is and js are the operations that need to be computed. See e.g. `_pdist_operations_generator_factory`
    :param n_operations: length of _operations_generator (as generators should not have __len__ method)
    :param n_processes: number of processes to use (defaults to maximum number of CPU cores)
    :param dtw_args: args to pass to dtw
    :param dtw_kwargs: kwargs to pass to dtw
    :return:
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

    # Create a shared memory array to store result
    # Do not lock it as the worker should make sure processes do not overlap the data
    result_buffer = Array(ctypes.c_double, n_operations, lock=False)

    # Create a shared memory buffer for the data array
    shape = three_dim_array.shape
    data_buffer = Array(ctypes.c_double, np.product(shape), lock=False)  # Allocate memory
    data_view = np.frombuffer(data_buffer).reshape(shape)  # Point numpy array to memory
    data_view[:] = three_dim_array  # Copy the contents into the new memory location

    number_of_slices = n_processes * 4

    buffer_size, remainder = divmod(n_operations, number_of_slices)

    scheduling_queue = Queue()

    # Split the data in slices
    for i in xrange(number_of_slices):
        start = i * buffer_size
        end = start + buffer_size
        scheduling_queue.put((start, end))
    if remainder:
        last_start = number_of_slices * buffer_size
        last_end = n_operations
        scheduling_queue.put((last_start, last_end))

    for i in xrange(n_processes):
        scheduling_queue.put(None)  # Add stop items to the queue so we know when its empty for sure

    # Create queue for exceptions
    exception_queue = Queue()

    processes = []
    for i in xrange(n_processes):
        p = Process(target=_parallel_dtw_worker,
                    args=(data_buffer, _operations_generator_factory, shape, result_buffer, scheduling_queue,
                          exception_queue, dtw_args, dtw_kwargs))
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
    n_items = three_dim_array.shape[0]
    number_of_combinations = combinations_count(n_items)

    return _parallel_dtw(three_dim_array, _pdist_operations_generator_factory, number_of_combinations,
                         n_processes=n_processes, *dtw_args, **dtw_kwargs)


def _path_calculation_worker(data_buffer, shape, prototypes_buffer, prototypes_shape,
                             scheduling_queue, results_queue, exception_queue,
                             dtw_args, dtw_kwargs):
    import sys
    import traceback
    import os

    try:
        pid = os.getpid()
        debug('PROCESS {0}: Spawned'.format(pid))
        data_view = np.frombuffer(data_buffer).reshape(shape)  # Point numpy array to memory

        prototypes_view = np.frombuffer(prototypes_buffer).reshape(prototypes_shape)
        answers = []
        while True:
            work_id = scheduling_queue.get()
            if work_id is None:
                break
            data_i, base_i = work_id
            x = data_view[data_i]
            base = prototypes_view[base_i]

            _, _, path = dtw_std(x, base, dist_only=False, *dtw_args, **dtw_kwargs)

            results_queue.put((work_id, path))

        debug('PROCESS {0}: finished computation'.format(pid))
        debug('PROCESS {0}: finished'.format(pid))

    except Exception as e:
        exception_queue.put(e)
        traceback.print_exc()
    except:
        e = sys.exc_info()[0]
        exception_queue.put(e)
        traceback.print_exc()

def parallel_dtw_paths(full_data, nodes, n_processes=None, *dtw_args, **dtw_kwargs):
    def _read_results_of_queue_till_empty(queue, ans_dict):
        debug('Queue full, reading results till can proceed further')
        number_of_answers = 0
        while number_of_answers < 100:
            try:
                (i, j), ps = queue.get(block=False)
            except Empty:
                debug('Answers queue empty, continuing')
                break
            ans_dict[node_id_lookup[j]][data_index[i]] = ps
            number_of_answers += 1

        debug('Got {0} answers back'.format(number_of_answers))
        return number_of_answers


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

    max_prototype_len = 0
    for node in nodes:
        prototype_len = len(node.prototype)
        if prototype_len > max_prototype_len:
            max_prototype_len = prototype_len
    ndims = full_data.values.shape[2]
    nan = [np.nan] * ndims
    prototypes = np.empty((len(nodes), max_prototype_len, ndims))
    for i, node in enumerate(nodes):
        prototype = node.prototype
        prototype_len = len(prototype)
        if prototype_len < max_prototype_len:
            padding = np.asarray([nan] * (max_prototype_len - prototype_len))
            padded_prototype = np.concatenate((prototype, padding))
            prototypes[i] = padded_prototype
        else:
            prototypes[i] = prototype

    # Create shared memory buffer fro prototypes array
    prototypes_shape = prototypes.shape
    prototypes_buffer = Array(ctypes.c_double, np.product(prototypes_shape), lock=False)
    prototypes_view = np.frombuffer(prototypes_buffer).reshape(prototypes_shape)
    prototypes_view[:] = prototypes

    # Generate lookup for indices
    data_index = full_data.items
    ix_lookup = {}
    for i, ix in enumerate(data_index):
        ix_lookup[ix] = i

    node_id_lookup = {}
    for j, node in enumerate(nodes):
        node_id_lookup[j] = node.id

    full_data = np.asarray(full_data)

    # Create a shared memory buffer for the data array
    shape = full_data.shape
    data_buffer = Array(ctypes.c_double, np.product(shape), lock=False)  # Allocate memory
    data_view = np.frombuffer(data_buffer).reshape(shape)  # Point numpy array to memory
    data_view[:] = full_data  # Copy the contents into the new memory location


    answers_queue = Queue()
    exception_queue = Queue()

    scheduling_queue = Queue()
    processes = []
    for i in xrange(n_processes):
        p = Process(target=_path_calculation_worker,
                    args=(data_buffer, shape, prototypes_buffer, prototypes_shape, scheduling_queue, answers_queue,
                          exception_queue, dtw_args, dtw_kwargs))
        processes.append(p)

    # Buffer to store results
    paths = {}
    for node in nodes:
        paths[node.id] = {}

    # Start all processes
    for p in processes:
        debug('Starting process')
        p.start()

    # Keep pushing things for the processes to work on
    n_items = 0

    number_of_answers = 0

    for j, node in enumerate(nodes):
        indices = node.index

        for ix in indices:
            while True:
                try:
                    scheduling_queue.put((ix_lookup[ix], j), block=False)
                    n_items += 1
                    break
                except Full:
                    # Whoops, our scheduling queue is full
                    debug('Full at {0} sent'.format(n_items))
                    number_of_answers += _read_results_of_queue_till_empty(answers_queue, paths)
                    continue

    # Stop all processes
    for i in xrange(n_processes):
        while True:
            try:
                scheduling_queue.put(None, block=False)
                break
            except Full:
                debug('Full while sending the Nones')
                number_of_answers += _read_results_of_queue_till_empty(answers_queue, paths)
                continue

    debug('Finished scheduling tasks, reading off results')
    while number_of_answers < n_items:
        (i, j), ps = answers_queue.get()
        paths[node_id_lookup[j]][data_index[i]] = ps
        number_of_answers += 1

    # Join all processes
    for p in processes:
        debug('Joining process')
        p.join()

    debug('Checking exceptions')
    try:
        exception = exception_queue.get(block=False)
        raise exception
    except Empty:
        pass

    return paths