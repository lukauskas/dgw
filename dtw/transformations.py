__author__ = 'saulius'
from distance import dtw_std
import pandas as pd
import numpy as np
from math import ceil, floor

def uniform_scaling_by_a_factor(sequence, scaling_factor):
    """
    Uses uniform scaling to scale the sequence by the scale_factor.
    .. [#yankov2007] D Yankov, E Keogh, J Medina, and B Chiu, "Detecting time series motifs under uniform scaling", 2007
    :param sequence:
    :param scaling_factor:
    :return:
    """
    scaling_factor = float(scaling_factor)
    current_len = len(sequence)
    rescaled_len = int(ceil(current_len / scaling_factor))

    # Using floor (unlike Yankov) as we are 0 based here
    rescaled_sequence = [sequence[int(floor(i*scaling_factor))] for i in range(rescaled_len)]

    return rescaled_sequence

def uniform_scaling_to_length(sequence, desired_length):
    current_len = len(sequence)
    scaling_factor = float(current_len) / desired_length

    rescaled_sequence = [sequence[int(floor(i*scaling_factor))] for i in range(desired_length)]
    return rescaled_sequence

def shrink_to_length(sequence, desired_length):
    current_length = len(sequence)

    # This is essentially how many points in the current sequence will be mapped to a single point in the newone
    shrink_factor = float(current_length) / desired_length

    new_sequence = np.empty(desired_length)
    for i in range(desired_length):
        start = int(floor(i * shrink_factor))
        end = int(floor((i+1) * shrink_factor))

        subset = sequence[start:end]
        new_sequence[i] = np.mean(subset)
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

def dtw_path_averaging(sequence_a, sequence_b, *args, **kwargs):
    """
    Averages the path computed between the two DTW sequences.
    Computes the DTW distance in order to do so.

    :param sequence_a: first sequence to be averaged
    :param sequence_b: second sequence to be averaged
    :param args: positional arguments passed into `dtw_std`
    :param kwargs: keyword arguments to be passed into  `dtw_std`
    :return:
    """

    distance, cost, path = dtw_std(sequence_a, sequence_b, dist_only=False, *args, **kwargs)

    path_base, path_other = path

    avg = np.array([(sequence_a[i] + sequence_b[j]) / 2.0 for i, j in zip(path_base, path_other)])

    return avg