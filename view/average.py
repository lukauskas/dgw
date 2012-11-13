'''
Created on 13 Nov 2012

@author: saulius
'''
from __future__ import print_function
from collections import defaultdict
from matplotlib import pyplot


def plot(data):
    
    sums   = defaultdict(lambda: 0)
    counts = defaultdict(lambda: 0)
    
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
    
    pyplot.plot(x,y)
    return (x,y)
  
    
    