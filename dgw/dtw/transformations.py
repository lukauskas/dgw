from itertools import izip

__author__ = 'saulius'
from dgw.dtw.distance import dtw_std
import pandas as pd
import numpy as np
from math import ceil, floor
from dgw.dtw.distance import _strip_nans

def uniform_scaling_to_length(sequence, desired_length):
    """
    Uniform scaling procedure, similar to the one provided in [#yankov2007]
    .. [#yankov2007] D Yankov, E Keogh, J Medina, and B Chiu, "Detecting time series motifs under uniform scaling", 2007
    :param sequence:
    :param desired_length:
    :return:
    """
    sequence = _strip_nans(sequence)
    current_len = len(sequence)
    if current_len == 0:
        raise ValueError('Empty sequence cannot be extended')
    elif desired_length == current_len:
        return sequence
    elif desired_length < current_len:
        raise ValueError('Desired length is smaller than current length: {0} < {1}'.format(desired_length, current_len))

    scaling_factor = float(current_len) / desired_length

    rescaled_sequence = [sequence[int(floor(i*scaling_factor))] for i in range(desired_length)]
    return rescaled_sequence

def uniform_shrinking_to_length(sequence, desired_length):
    sequence = _strip_nans(sequence)
    sequence = np.asarray(sequence, dtype=float)

    current_length = len(sequence)

    if current_length == 0:
        raise ValueError('Cannot shrink sequence of length 0')
    elif current_length < desired_length:
        raise ValueError('Desired length greater than current length: {0} > {1}'.format(desired_length, current_length))
    elif current_length == desired_length:
        return sequence

    if desired_length <= 0:
        raise ValueError('Invalid length desired: {0}'.format(desired_length))

    # This is essentially how many points in the current sequence will be mapped to a single point in the newone
    shrink_factor = float(current_length) / desired_length

    try:
        ndim = sequence.shape[1]
    except IndexError:
        ndim = 0
    if ndim == 0:
        new_sequence = np.empty(desired_length)
    else:
        new_sequence = np.empty((desired_length, ndim))

    for i in range(desired_length):
        start = i * shrink_factor
        end = (i + 1) * shrink_factor

        s = 0
        d = 0

        left_bound = int(floor(start))

        if start == left_bound:
            left_bound_input = 1
        else:
            left_bound_input = ceil(start) - start

        if left_bound_input > 0:
            s += sequence[left_bound] * left_bound_input
            d += left_bound_input

        right_bound = int(floor(end))
        right_bound_input = end - floor(end)

        if right_bound_input > 0:
            print right_bound, right_bound_input, sequence[right_bound], sequence[right_bound] * right_bound_input
            s += sequence[right_bound] * right_bound_input
            d += right_bound_input

        for j in xrange(left_bound + 1, right_bound):
            s += sequence[j]
            d += 1.0

        assert(abs(d - shrink_factor) < 0.000001)

        new_sequence[i] = s / d

    return new_sequence




def dtw_projection(sequence, base_sequence, *args, **kwargs):
    """
    Projects given sequence onto a base time series using Dynamic Time Warping
    :param sequence: the sequence that will be projected onto base_sequence
    :param base_sequence: base sequence to project onto
    :return: new time series of length base containing x projected on it
    """
    distance, cost, path = dtw_std(base_sequence, sequence, dist_only=False, *args, **kwargs)

    path_base, path_other = path

    current_sums = np.zeros(len(base_sequence))
    current_counts = np.zeros(len(base_sequence))

    current_sums[np.isnan(base_sequence)] = np.nan
    current_counts[np.isnan(base_sequence)] = np.nan

    for mapped_i, i in zip(path_base, path_other):
        # Go through the path and sum all points that map to the base location i together
        current_sums[mapped_i] += sequence[i]
        current_counts[mapped_i] += 1

    current_average = current_sums / current_counts

    return current_average

def dtw_projection_multi(sequences, base, *args, **kwargs):
    """
    Projects given sequences onto the base.

    :param base:
    :type base: array-like
    :param data_frame:
    :type data_frame: pd.DataFrame or np.ndarray
    :param args:
    :param kwargs:
    :return:
    :rtype: pd.DataFrame or np.ndarray
    """

    def _projection(sequence):
        return dtw_projection(sequence, base, *args, **kwargs)

    index = None
    if isinstance(sequences, pd.DataFrame):
        index = sequences.index
        sequences = sequences.values

    projected_sequences = map(_projection, sequences)

    if index is not None:
        return pd.DataFrame(projected_sequences, index=index)
    else:
        return projected_sequences

def dtw_path_averaging(sequence_a, sequence_b, path=None, *args, **kwargs):
    """
    Averages the path computed between the two DTW sequences.
    Computes the DTW distance in order to do so.

    :param sequence_a: first sequence to be averaged
    :param sequence_b: second sequence to be averaged
    :param path: (optional) computed mapped path between the sequences. Will be computed if not provided
    :param args: positional arguments passed into `dtw_std`
    :param kwargs: keyword arguments to be passed into  `dtw_std`
    :return:
    """
    sequence_a = np.asarray(sequence_a, dtype=float)
    sequence_b = np.asarray(sequence_b, dtype=float)

    if path is None:
        distance, cost, path = dtw_std(sequence_a, sequence_b, dist_only=False, *args, **kwargs)

    path_base, path_other = path

    avg = np.array([(sequence_a[i] + sequence_b[j]) / 2.0 for i, j in zip(path_base, path_other)])

    return avg

def sdtw_averaging(sequence_a, sequence_b, weight_a, weight_b, path=None, *args, **kwargs):
    """
    Implements Scaled Dynamic Time Warping Path Averaging as described in [#Niennattrakul:2009ep]

    .. [#Niennattrakul:2009ep] Vit Niennattrakul and Chotirat Ann Ratanamahatana "Shape averaging under Time Warping",
       2009 6th International Conference on Electrical Engineering/Electronics, Computer, Telecommunications and Information Technology (ECTI-CON)
    :param sequence_a: sequence A
    :param sequence_b: sequence B
    :param weight_a: weight of sequence A
    :param weight_b: weight of sequence B
    :param path: computed DTW path. Will be calculated automatically if not provided
    :param args: args passed to `dtw_std`
    :param kwargs: kwargs passed to `dtw_std`
    :return:
    """
    sequence_a = np.asarray(sequence_a, dtype=float)
    sequence_b = np.asarray(sequence_b, dtype=float)

    if path is None:
        distance, cost, path = dtw_std(sequence_a, sequence_b, dist_only=False, *args, **kwargs)

    path = izip(path[0], path[1])  # Rezip this for easier traversal

    averaged_path = []

    prev = None

    diagonal_coefficient = int((weight_a + weight_b) / 2.0)  # The paper does not explicitly say how to round this
    for a,b in path:

        item = float(weight_a * sequence_a[a] + weight_b * sequence_b[b]) / (weight_a + weight_b)
        if prev is None:
            extension_coefficient = diagonal_coefficient
        else:
            # TODO: remove assertions
            if prev[0] == a:  # The path moved from (i,j-1) to (i,j)
                assert(prev[1] + 1 == b)
                extension_coefficient = weight_a
            elif prev[1] == b:  # The path moved from (i-1,j) to (i,j)
                assert(prev[0] + 1 == a)
                extension_coefficient = weight_b
            else:  # Path moved diagonally from (i-1,j-1) to (i,j)
                assert(prev[0] + 1 == a)
                assert(prev[1] + 1 == b)
                extension_coefficient = diagonal_coefficient

        new_items = [item] * extension_coefficient
        averaged_path.extend(new_items)
        prev = (a, b)

    return uniform_shrinking_to_length(averaged_path, max(len(_strip_nans(sequence_a)), len(_strip_nans(sequence_b))))















