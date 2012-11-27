'''
Created on 13 Nov 2012

@author: saulius
'''
import pysam
from peak import Peak
from itertools import izip
from tss_peak_analyser import read_tss, KNOWNGENES_FILE
from view.average import plot as plot_average
from matplotlib import pyplot
from cluster.distance.dtw import dtw_distance_matrix
from cluster.cluster import hierarchical_cluster
import math
from cache import cached
import os
from mlpy import dtw_std
from view.heatmap import plot as plot_heatmap

from math import floor, ceil
import numpy


def parse_peak_from_bam_line(line, resolution = 1):
    '''
    @rtype: Peak
    '''
    line = line.strip("\n")
    line = line.split("\t")
    
    chromosome = line[0]
    start      = resolution * int(floor(float(line[1]) / resolution))
    end        = resolution * int(ceil(float(line[2]) / resolution))
    
    return Peak(chromosome, start, end, resolution=resolution)

def get_peak_data(samfile, peak):
    '''
    @param samfile: samfile to query
    @type samfile: pysam.Samfile
    @param peak: peak to get the data for
    @type peak: Peak
    @rtype: list
    '''
    
    data_arr = []
    
    # Note this has to be peak.end - 1 
    # because samfile.pileup treats ends inclusively
    data_len = (peak.end - peak.start) / peak.resolution
#    return numpy.random.rand(data_len)
    
    for i in xrange(data_len):
        start = peak.start + i * peak.resolution
        end   = start + peak.resolution - 1 # samfile treats ends inclusively
        assert(end > start)
        
        count = samfile.count(peak.chromosome, start, end)
        data_arr.append(count)
    
    assert(len(data_arr) > 0)
    assert(len(data_arr) == data_len)
    return data_arr    

@cached
def get_peaks(peaks_file, alignments_file, tss_set, resolution=1):
    
    samfile = pysam.Samfile(alignments_file, "rb")
    peaks_handle = open(peaks_file, "r")
    
    counter = 0
    peaks = []
    for line in peaks_handle:

        peak = parse_peak_from_bam_line(line, resolution)
        
        peak.find_interesting_points_from_set(tss_set)
        if len(peak.points_of_interest) != 1:
            continue
        peak.data = get_peak_data(samfile, peak)
        
        peaks.append(peak)
        
        counter += 1
        if counter > 100:
            break

    samfile.close()
    peaks_handle.close()
    
    return peaks

if __name__ == '__main__':
    tss_set = read_tss(KNOWNGENES_FILE)
    
    import sys
    
    print("Obtaining data")
    peaks = list(get_peaks(*sys.argv[1:], tss_set=tss_set, resolution=25))
    print len(peaks)
    
    data_extract_func = lambda peak: peak.data_relative_to_start()
    data = map(lambda x : map(lambda y : y[1], sorted(data_extract_func(x))), peaks)
    N = len(data)
   
    
    print("Clustering")
    THRESHOLD = 10
    distance_func = cached(dtw_distance_matrix, os.path.basename(sys.argv[1]))
    cluster_assignments = hierarchical_cluster(data, distance_func, THRESHOLD)
    print("Number of clusters")
    cluster_assignments = list(cluster_assignments)
    number_of_clusters = max(cluster_assignments)
    
    ROWS = int(math.sqrt(number_of_clusters))
    COLS = ROWS
    while (ROWS*COLS < number_of_clusters):
        ROWS += 1
        
    print('Rows: {0}, cols: {1}'.format(ROWS,COLS))
    
    pyplot.figure(1)
    for i in range(1, number_of_clusters+1):
        pyplot.subplot(ROWS, COLS, i)
        cluster_data = [ p for c,p in izip(cluster_assignments, peaks) if c == i]
        print("Cluster {0}: {1}".format(i, len(cluster_data)))
        plot_average(cluster_data, data_extract_func)
    
    pyplot.savefig('averages.png')
    
    fig = pyplot.figure(2)
    for i in range(1, number_of_clusters+1):
        cluster_data = [ p for c,p in izip(cluster_assignments, peaks) if c == i]
        pyplot.subplot(ROWS,COLS,i)
        plot_heatmap(cluster_data, data_extract_func)
    
    pyplot.savefig('heatmaps.png')
   
#    fig = pyplot.figure()
#    ax = fig.add_subplot(1, 2, 1)
#    data_extract_func = lambda peak: peak.data_relative_to_start()
#    plot_average(peaks, data_extract_func)
#    ax = fig.add_subplot(1,2,2)
#    plot_heatmap(peaks, data_extract_func, ax)
#    pyplot.savefig("out.png")
#    pyplot.show()

#        
    #print("Number of clusters: {0}".format(max(cluster_assignments)))
    
     
      