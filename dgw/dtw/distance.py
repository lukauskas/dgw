from dgw._mlpy.dtw import dtw_std as mlpy_dtw_std
import numpy as np
from dgw.dtw.scaling import uniform_scaling_to_length
from dgw.dtw.utilities import _strip_nans, no_nans_len, reverse_sequence


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
            scale_first=False, *args, **kwargs):
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
    :param scale_first: If set to true, the shorte sequence will be scaled to the length of the longer sequence before DTW
    :param kwargs:
    :return:
    """
    def _normalise(ans, max_len):
        if normalise:
            return ans / max_len
        else:
            return ans

    def _scaled_path(path, scaling_path, flip_paths):
        path_x = np.asarray([scaling_path[i] for i in path[0]])
        path_y = path[1]

        if flip_paths:
            path = (path_y, path_x)
        else:
            path = (path_x, path_y)

        return path

    def _reverse_path(path):
        n = path.max()
        path = n - path
        return path


    x = np.asarray(x, dtype=np.float)
    y = np.asarray(y, dtype=np.float)

    x = _strip_nans(x)
    y = _strip_nans(y)

    max_len = max(len(x), len(y))
    if scale_first:
        if len(x) >= len(y):
            x, y = y, x
            flip_paths = True
        else:
            flip_paths = False

        x, scaling_path = uniform_scaling_to_length(x, len(y), output_scaling_path=True)

    regular_ans = mlpy_dtw_std(x, y, metric=metric, dist_only=dist_only, constraint=constraint, k=k, *args, **kwargs)
    if not try_reverse:
        if dist_only:
            return _normalise(regular_ans, max_len)
        else:
            dist, cost, path = regular_ans
            dist = _normalise(dist, max_len)

            if scale_first:
                path = _scaled_path(path, scaling_path, flip_paths)

            return dist, cost, path
    else:
        reverse_ans = mlpy_dtw_std(reverse_sequence(x), y, metric=metric, dist_only=dist_only, constraint=constraint, k=k, *args, **kwargs)
        if dist_only:
            return _normalise(min(regular_ans, reverse_ans), max_len)
        elif reverse_ans[0] >= regular_ans[0]:
            dist, cost, path = regular_ans
            if scale_first:
                path = _scaled_path(path, scaling_path, flip_paths)
            return _normalise(dist, max_len), cost, path
        else:  # dist_only = False and reverse_ans is smaller
            dist, cost, path = reverse_ans
            path_rev = (_reverse_path(path[0]), path[1])

            if scale_first:
                path_rev = _scaled_path(path_rev, scaling_path, flip_paths)

            cost = np.fliplr(cost)
            return _normalise(dist, max_len), cost, path_rev

def dtw_path_is_reversed(path):
    """
    Returns true if DTW path is reversed

    :param path:
    :return:
    """
    # Just need to check whether first point in the first sequence path is zero.
    return path[0][0] != 0

def warping_conservation_vector(warping_path):
    """
    Computes warping conservation vector for the warping path given.
    This vector is always of the length n-1 where n is the length of the second sequence (usually the base sequence).
    This new vector contains 1 for whenever the base sequence is conserved (the path moves diagonally),
    and contains 0 otherwise.

    :param warping_path:
    :return:
    """
    path_a, path_b = warping_path

    n = max(path_b[0], path_b[-1]) + 1  # Number of points on the second sequence is, either the last point
                                           # .. or the first point (if reversed), plus one

    conservation_vector = np.zeros(n-1)

    prev_i = None
    prev_j = None
    diag_path_taken_until_current_point = True
    for i, j in zip(path_a, path_b):
        if prev_j is None:
            prev_j = j
            prev_i = i
            diag_path_taken_until_current_point = True
            continue
        elif j == prev_j:
            diag_path_taken_until_current_point = False
            prev_i = i
        else:
            if diag_path_taken_until_current_point and prev_i != i:
                conservation_vector[min(prev_j, j)] = 1

            diag_path_taken_until_current_point = True
            prev_j = j
            prev_i = i

    first_one_in_set = -1
    number_of_ones_in_set = 0
    for i in range(n-1):
        if conservation_vector[i] == 0 and first_one_in_set > -1:
            conservation_vector[first_one_in_set:i] = number_of_ones_in_set
            number_of_ones_in_set = 0
            first_one_in_set = -1
        elif conservation_vector[i] == 1:
            if number_of_ones_in_set == 0:
                first_one_in_set = i
            number_of_ones_in_set += 1

    if number_of_ones_in_set > 0 and first_one_in_set > -1:
        conservation_vector[first_one_in_set:] = number_of_ones_in_set

    return conservation_vector
