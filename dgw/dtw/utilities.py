import numpy as np

__author__ = 'saulius'


def _strip_nans(sequence):
    '''
        Strips NaNs that are padded to the right of each peak if they are of unequal length
    :param sequence:
    :return:
    '''
    sequence = np.asarray(sequence)

    try:
        lookup = np.all(np.isnan(sequence), axis=1)
    except ValueError:
        # Will get thrown for one-dimensional arrays
        return sequence[~np.isnan(sequence)]

    sequence = sequence[~lookup]

    if np.any(np.isnan(sequence)):
        raise ValueError('Inconsistent NaNs between dimensions')

    return sequence


def no_nans_len(sequence):
    """
    Returns length of the sequence after removing all the nans from it.

    :param sequence:
    :return:
    """
    return len(_strip_nans(sequence))

def reverse_sequence(sequence):
    no_nans_sequence = _strip_nans(sequence)
    rev = no_nans_sequence[::-1]

    padding = len(sequence) - len(no_nans_sequence)

    try:
        ndim = no_nans_sequence.shape[1]
    except IndexError:
        return np.concatenate((rev, [np.nan] * padding))

    if padding:
        nans = [[np.nan] * ndim] * padding
        return np.concatenate((rev, nans))
    else:
        return rev
