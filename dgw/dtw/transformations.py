from itertools import izip
from logging import debug
import pandas as pd
import numpy as np
from math import ceil, floor, fabs
from dgw.data.containers import AlignmentsData

from distance import dtw_std, _strip_nans, no_nans_len

def points_mapped_to(point_on_original_sequence, dtw_path, sequence_a=True):
    """
    Returns all indices on the sequence b that are mapped to the index on sequence_a provided.

    :param point_on_original_sequence: index of a point on original sequence
    :param dtw_path: DTW warping path between `sequence_a` and `sequence_b`
    :param sequence_a: whether the query point is on `sequence_a` or not. (Will swap sequences otherwise)
    :return:
    """
    if sequence_a:
        path_ours, path_theirs = dtw_path
    else:
        path_theirs, path_ours = dtw_path

    path_indices = np.nonzero(path_ours == point_on_original_sequence)
    return path_theirs[path_indices]

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
    EPSILON = 1e-6

    sequence = np.asarray(sequence, dtype=float)
    sequence = _strip_nans(sequence)

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

        if fabs(start - left_bound) <= EPSILON:
            left_bound_input = 1
        else:
            left_bound_input = ceil(start) - start

        if left_bound_input >= EPSILON:
            s += sequence[left_bound] * left_bound_input
            d += left_bound_input

        right_bound = int(floor(end))
        right_bound_input = end - floor(end)

        if right_bound_input >= EPSILON:  # Epsilon to prevent rounding errors interfering
            s += sequence[right_bound] * right_bound_input
            d += right_bound_input

        for j in xrange(left_bound + 1, right_bound):
            s += sequence[j]
            d += 1.0

        assert(abs(d - shrink_factor) < 0.000001)

        new_sequence[i] = s / d

    return new_sequence




def dtw_projection(sequence, base_sequence, dtw_function=dtw_std, path=None):
    """
    Projects given sequence onto a base time series using Dynamic Time Warping
    :param sequence: the sequence that will be projected onto base_sequence
    :param base_sequence: base sequence to project onto
    :param dtw_function: DTW function to compute the path if it is set to None
    :param path: the pre-computed DTW warping path between sequence and base sequence.
    :return: new time series of length base containing x projected on it
    """
    base_sequence = np.asarray(base_sequence)
    sequence = np.asarray(sequence)

    if not path:
        distance, cost, path = dtw_function(sequence, base_sequence, dist_only=False)

    path_other, path_base = path

    current_sums = np.zeros(base_sequence.shape)
    current_counts = np.zeros(base_sequence.shape)

    nnl = no_nans_len(base_sequence)

    try:
        filler = [np.nan] * base_sequence.shape[1]
    except IndexError:
        filler = np.nan

    for i in range(nnl, len(base_sequence)):
        current_sums[i] = filler
        current_counts[i] = filler

    for mapped_i, i in zip(path_base, path_other):
        # Go through the path and sum all points that map to the base location i together
        current_sums[mapped_i] += sequence[i]
        current_counts[mapped_i] += 1

    current_average = current_sums / current_counts

    # Append NaNs as needed
    nans_count = len(base_sequence) - len(base_sequence)
    if nans_count:
        current_average = np.concatenate((current_average,
                                          [[np.nan] * base_sequence.shape[-1]] * (nans_count)))

    return current_average

def dtw_projection_multi(alignments, base, *args, **kwargs):
    """
    DTW Projects given alignments onto the given base
    """
    def _projection(sequence):
        return dtw_projection(sequence, base, *args, **kwargs)

    index = None
    if isinstance(alignments, pd.DataFrame):
        index = alignments.index
        sequences = alignments.values

    new_data = {}

    for index in alignments.items:
        item = alignments.ix[index]
        projection = dtw_projection(item, base, *args, **kwargs)
        new_data[index] = pd.DataFrame(projection, index=range(len(base)), columns=item.columns)

    return AlignmentsData(pd.Panel(new_data))

def dtw_path_averaging(sequence_a, sequence_b, weight_a=1, weight_b=1, path=None, shrink=True, dtw_function=dtw_std):
    """
    Averages the path computed between the two DTW sequences.
    Computes the DTW distance in order to do so.

    :param sequence_a: first sequence to be averaged
    :param sequence_b: second sequence to be averaged
    :param weight_a: weight of first sequence
    :param weight_b: weight of the second sequence
    :param path: computed mapped path between the sequences. Will be computed using `dtw_function` if not provided
    :param shrink: if set to treu the data will be shrinked to length of maximum sequence
    :param dtw_function: function that computes DTW path between two sequences, e.g. see `parametrised_dtw_wrapper`
    :return:
    """
    sequence_a = np.asarray(sequence_a, dtype=float)
    sequence_b = np.asarray(sequence_b, dtype=float)

    if path is None:
        distance, cost, path = dtw_function(sequence_a, sequence_b, dist_only=False)

    path_base, path_other = path

    averaged_path = np.array([(sequence_a[i] * weight_a + sequence_b[j] * weight_b) / (weight_a + weight_b) for i, j in zip(path_base, path_other)])

    if shrink:
        averaged_path = uniform_shrinking_to_length(averaged_path, max(len(_strip_nans(sequence_a)),
                                                                       len(_strip_nans(sequence_b))))

    return averaged_path

def sdtw_averaging(sequence_a, sequence_b, weight_a, weight_b, path=None, shrink=True, dtw_function=dtw_std):
    """
    Implements Scaled Dynamic Time Warping Path Averaging as described in [#Niennattrakul:2009ep]

    .. [#Niennattrakul:2009ep] Vit Niennattrakul and Chotirat Ann Ratanamahatana "Shape averaging under Time Warping",
       2009 6th International Conference on Electrical Engineering/Electronics, Computer, Telecommunications and Information Technology (ECTI-CON)

    :param sequence_a: sequence A
    :param sequence_b: sequence B
    :param weight_a: weight of sequence A
    :param weight_b: weight of sequence B
    :param path: computed mapped path between the sequences. Will be computed using `dtw_function` if not provided
    :param shrink: if set to true the data will be shrinked to the length of maximum seq
    :return:
    """
    sequence_a = np.asarray(sequence_a, dtype=float)
    sequence_b = np.asarray(sequence_b, dtype=float)

    if path is None:
        distance, cost, path = dtw_function(sequence_a, sequence_b, dist_only=False)

    path = izip(path[0], path[1])  # Rezip this for easier traversal

    averaged_path = []

    prev = None

    diagonal_coefficient = int((weight_a + weight_b) / 2.0)  # The paper does not explicitly say how to round this
    for a,b in path:

        item = (weight_a * sequence_a[a] + weight_b * sequence_b[b]) / (weight_a + weight_b)
        if prev is None:
            extension_coefficient = diagonal_coefficient
        else:
            if prev[0] == a:  # The path moved from (i,j-1) to (i,j)
                # assert(prev[1] + 1 == b)
                extension_coefficient = weight_a
            elif prev[1] == b:  # The path moved from (i-1,j) to (i,j)
                # assert(prev[0] + 1 == a)
                extension_coefficient = weight_b
            else:  # Path moved diagonally from (i-1,j-1) to (i,j)
                # assert(prev[0] + 1 == a)
                # assert(prev[1] + 1 == b)
                extension_coefficient = diagonal_coefficient

        new_items = [item] * extension_coefficient
        averaged_path.extend(new_items)
        prev = (a, b)

    averaged_path = np.asarray(averaged_path, dtype=float)
    if shrink:
        averaged_path = uniform_shrinking_to_length(averaged_path, max(len(_strip_nans(sequence_a)),
                                                                       len(_strip_nans(sequence_b))))
    return averaged_path














