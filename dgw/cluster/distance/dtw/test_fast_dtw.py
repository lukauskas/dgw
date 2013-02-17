import random

import numpy
from scipy.spatial.distance import cityblock

from dgw.dtw import fast_dtw
from dgw import dtw


def main(count, seed = None):
    if seed is not None:
        numpy.random.seed(seed)
        random.seed(seed)
    
    MIN_LEN = 5
    MAX_LEN = 40
    test_data = []
    for _ in range(count):
        pair = (numpy.random.rand(random.randint(MIN_LEN, MAX_LEN)), numpy.random.rand(random.randint(MIN_LEN,MAX_LEN)))
        
        test_data.append(pair)
    
    for a, b in test_data:
        print '{0!r}, {1!r}:'.format([ x for x in a], [x for x in b])
        dtw_std,_ = dtw(a,b, cityblock)
        dtw_fast,_ = fast_dtw(a,b, cityblock)
        
        ratio = float(dtw_fast) / dtw_std
        
        print '>> {ratio} ({dtw_fast}/{dtw_std})'.format(ratio=ratio,
                                                                 dtw_fast=dtw_fast,
                                                                 dtw_std=dtw_std)

if __name__ == '__main__':
    main(1000, 10)
    
    