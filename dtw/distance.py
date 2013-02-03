__author__ = 'saulius'

from mlpy.dtw import dtw_std as mlpy_dtw_std
import numpy as np

def _strip_nans(sequence):
    '''
        Strips NaNs that are padded to the right of each peak if they are of unequal length
    :param sequence:
    :return:
    '''
    return sequence[~np.isnan(sequence)]

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

    if metric == 'sqeuclidean':
        squared = True
    elif metric == 'euclidean':
        squared = False
    else:
        raise ValueError('Unsupported metric provided: {0!r}'.format(metric))

    return mlpy_dtw_std(x, y, squared=squared, *args, **kwargs)


