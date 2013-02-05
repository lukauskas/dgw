__author__ = 'saulius'
from distance import dtw_std
import pandas as pd
import numpy as np

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