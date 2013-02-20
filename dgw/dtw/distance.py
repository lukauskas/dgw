from logging import debug

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

def dtw_std(x, y, metric='sqeuclidean', dist_only=True, constraint=None, k=None, try_reverse=True, *args, **kwargs):
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

    regular_ans = mlpy_dtw_std(x, y, metric=metric, dist_only=dist_only, constraint=constraint, k=k, *args, **kwargs)
    if not try_reverse:
        return regular_ans
    else:
        reverse_ans = mlpy_dtw_std(x[::-1], y, metric=metric, dist_only=dist_only, constraint=constraint, k=k, *args, **kwargs)
        if dist_only:
            return min(regular_ans, reverse_ans)
        elif reverse_ans[0] >= regular_ans[0]:
            return regular_ans
        else:  # dist_only = False and reverse_ans is smaller
            dist, cost, path = reverse_ans
            path_rev = (path[0][::-1], path[1])

            # TODO: Reverse cost matrix here
            cost = None
            return dist, cost, path_rev
