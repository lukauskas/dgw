from math import floor, fabs, ceil
import numpy as np
from dgw.dtw.utilities import _strip_nans


def uniform_scaling_to_length(sequence, desired_length, output_scaling_path=False):
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
        if output_scaling_path:
            return sequence, range(desired_length)
        else:
            return sequence
    elif desired_length < current_len:
        raise ValueError('Desired length is smaller than current length: {0} < {1}'.format(desired_length, current_len))

    scaling_factor = float(current_len) / desired_length

    rescaled_sequence = [sequence[int(floor(i*scaling_factor))] for i in range(desired_length)]

    if output_scaling_path:
        return rescaled_sequence, np.asarray([int(floor(i * scaling_factor)) for i in range(desired_length)])
    else:
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