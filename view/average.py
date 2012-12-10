'''
Created on 13 Nov 2012

@author: saulius
'''
from __future__ import print_function
from collections import defaultdict
from matplotlib import pyplot


def plot(peaks):

    sums   = defaultdict(lambda: 0)
    counts = defaultdict(lambda: 0)
    data = map(data_extract_func, peaks)
    
    for peak_data in data:
        for pos, n in peak_data:
            sums[pos]   += n
            counts[pos] += 1
        
    x = []
    y = []
    
    for i in sorted(sums):
        x.append(i)
        n = sums[i]
        count = counts[i]
        if count == 0:
            print(count)
        y.append(float(n) / count)
    
    pyplot.plot(x,y, figure=figure)
    return (x,y)
  
    
    