from logging import debug
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

def parametrised_dtw_wrapper(*dtw_args, **dtw_kwargs):
    """
    Returns a wrapper around DTW function with args dtw_args and kwargs dtw_kwargs
    :param dtw_args: positional parameters of DTW
    :param dtw_kwargs: keyword parameters of DTW
    :return:
    """

    def f(x, y, dist_only=False):
        return dtw_std(x, y, dist_only=dist_only, *dtw_args, **dtw_kwargs)

    return f

def dtw_std(x, y, metric='sqeuclidean', dist_only=True, constraint=None, k=None, try_reverse=True, normalise=False,
            *args, **kwargs):
    """
    Wrapper arround MLPY's dtw_std that supports cleaning up of NaNs, and reversing of strings.
    :param x:
    :param y:
    :param metric: dtw metric to use `sqeuclidean`, `euclidean` or `cosine`
    :param dist_only: return distance only
    :param constraint: constraint of dtw (try `None` or `'slanted_band'`
    :param k: parameter k needed for slanted band constraint
    :param try_reverse: Will try reversing one sequence as to get a better distance
    :param normalise: If set to true, distance will be divided from the length of the longer sequence
    :param kwargs:
    :return:
    """
    def _normalise(ans):
        if normalise:
            return ans / max(no_nans_len(x), no_nans_len(y))
        else:
            return ans

    x = np.asarray(x, dtype=np.float)
    y = np.asarray(y, dtype=np.float)

    x = _strip_nans(x)
    y = _strip_nans(y)

    regular_ans = mlpy_dtw_std(x, y, metric=metric, dist_only=dist_only, constraint=constraint, k=k, *args, **kwargs)
    if not try_reverse:
        return _normalise(regular_ans)
    else:
        reverse_ans = mlpy_dtw_std(x[::-1], y, metric=metric, dist_only=dist_only, constraint=constraint, k=k, *args, **kwargs)
        if dist_only:
            return _normalise(min(regular_ans, reverse_ans))
        elif reverse_ans[0] >= regular_ans[0]:
            dist, cost, path = regular_ans
            return _normalise(dist), cost, path
        else:  # dist_only = False and reverse_ans is smaller
            dist, cost, path = reverse_ans
            path_rev = (path[0][::-1], path[1])

            # TODO: Reverse cost matrix here
            cost = None
            return _normalise(dist), cost, path_rev
