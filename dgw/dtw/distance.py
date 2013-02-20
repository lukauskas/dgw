__author__ = 'saulius'

from mlpy.dtw import dtw_std as mlpy_dtw_std
import numpy as np

def _strip_nans(sequence):
    '''
        Strips NaNs that are padded to the right of each peak if they are of unequal length
    :param sequence:
    :return:
    '''

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

def dtw_std(x, y, metric='sqeuclidean', *args, **kwargs):
    '''
        Wrapper around mlpy's dtw_std that first strips all NaNs out of the data.

    @param x:
    @param y:
    @param args:
    @param kwargs:
    @return:
    '''
    x = np.asarray(x, dtype=np.float)
    y = np.asarray(y, dtype=np.float)

    x = _strip_nans(x)
    y = _strip_nans(y)

    return mlpy_dtw_std(x, y, metric=metric, *args, **kwargs)

def dtw_with_reversing(x, y, dist_only=True, *args, **kwargs):

    # reverse of x
    rev_x = x[::-1]

    if dist_only:
        dist = dtw_std(x, y, dist_only=True, *args, **kwargs)
        dist_rev = dtw_std(rev_x, y, dist_only=True, *args, **kwargs)

        if dist >= dist_rev:
            return dist, False
        else:
            return dist_rev, True

    else:
        dist, cost, path = dtw_std(x, y, dist_only=False, *args, **kwargs)
        dist_rev, cost_rev, path_rev = dtw_std(rev_x, y, dist_only=True, *args, **kwargs)

        if dist >= dist_rev:
            return (dist, cost, path), False
        else:
            cost_rev = None  # TODO: reverse cost matrix here
            path_rev = (path_rev[0][::-1], path_rev[1][::-1])
            return (dist_rev, cost_rev, path_rev), True



