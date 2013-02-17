'''
Created on 13 Nov 2012

@author: saulius
'''
from itertools import izip
import math
import os
from math import floor, ceil
import random

import pysam
from matplotlib import pyplot

from peak import Peak
from tss_peak_analyser import read_tss, KNOWNGENES_FILE
from view.average import plot as plot_average
from cluster.distance.dtw import dtw_distance_matrix
from cluster.cluster import hierarchical_cluster
from cache import cached
from dgw.data.visualisation.heatmap import plot as plot_heatmap


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
        if random.randint(1, 20000) > 7000:
            continue
        if counter == 7000:
            break

        peak = parse_peak_from_bam_line(line, resolution)
        
        peak.find_interesting_points_from_set(tss_set)
#        if len(peak.points_of_interest) != 1:
#            continue
        peak.data = get_peak_data(samfile, peak)
        
        peaks.append(peak)
        
        counter += 1
#        if counter > 2000:
#            break

    samfile.close()
    peaks_handle.close()
    
    return peaks

@cached
def get_near_tss_regions(alignments_file, tss_set, offset=2000, resolution=1):
    
    samfile = pysam.Samfile(alignments_file, "rb")

    peaks = []
    while (2*offset + 1) % resolution != 0:
        offset += 1
    print "Offset used: {0} bp".format(offset)
    
    
    for (chrom, pos) in tss_set:
        start      = pos-offset
        end        = pos+offset+1 # we want to include the final BP
        
        peak = Peak(chrom, start, end, resolution=resolution)
        peak.add_point_of_interest(pos)
        try:
            peak.data = get_peak_data(samfile, peak)
        except ValueError:
            continue
        
        peaks.append(peak)
    
    samfile.close()
    return peaks
        

def plot_clusters(cluster_assignments, plotting_function):
    
    number_of_clusters = max(cluster_assignments)
    
    print("Number of clusters: {0}".format(number_of_clusters))
    ROWS = int(math.sqrt(number_of_clusters))
    COLS = ROWS
    while (ROWS*COLS < number_of_clusters):
        ROWS += 1
        
    pyplot.figure()
    for i in range(1, number_of_clusters+1):
        pyplot.subplot(ROWS, COLS, i)
        cluster_data = [ p for c,p in izip(cluster_assignments, peaks) if c == i]
        plotting_function(cluster_data)
 
def print_cluster_sizes(cluster_assignments):
    
    
    number_of_clusters = max(cluster_assignments)
    print 'Number of clusters: {0}'.format(number_of_clusters)
    
    for i in range(1, number_of_clusters+1):
        cluster_data = [ p for c,p in izip(cluster_assignments, peaks) if c == i]
        print 'Cluster #{0}: {1}'.format(i, len(cluster_data))
        
if __name__ == '__main__':
    tss_set = read_tss(KNOWNGENES_FILE)
    random.seed(42)
    
    import sys
    
    print("Obtaining data")
    peaks = list(get_peaks(*sys.argv[1:], tss_set=tss_set, resolution=25))
    #peaks = list(get_near_tss_regions(sys.argv[2], tss_set=tss_set, offset=2000, resolution=25))
    
    # Get a random choice of peaks
    #peaks = random.sample(peaks, 5000)
    #peaks = peaks[:1000]
    print len(peaks)
    
    data_extract_func = lambda peak: peak.data_relative_to_start()
    data = map(lambda x : map(lambda y : y[1], sorted(data_extract_func(x))), peaks)
    N = len(data)
   
    
    #data_extract_func2 = lambda peak: peak.data_relative_to_point(list(peak.points_of_interest)[0])
    data_extract_func2 = data_extract_func
    print("Clustering")
    THRESHOLD = 1.2
    distance_func = cached(dtw_distance_matrix, os.path.basename(sys.argv[1]))
    cluster_assignments, z = hierarchical_cluster(data, distance_func, THRESHOLD, criterion="inconsistent", depth=2)
    
    cluster_assignments = list(cluster_assignments)
    print_cluster_sizes(cluster_assignments)
    
    plot_clusters(cluster_assignments, lambda x : plot_average(x, data_extract_func2))
    pyplot.savefig('averages.png')
    plot_clusters(cluster_assignments, lambda x : plot_heatmap(x, data_extract_func2))
    pyplot.savefig('heatmaps.png')
   
#    fig = pyplot.figure()
#    ax = fig.add_subplot(1, 2, 1)
#    plot_average(peaks, data_extract_func)
#    ax = fig.add_subplot(1,2,2)
#    plot_heatmap(peaks, data_extract_func)
#    pyplot.savefig("out.png")
#    pyplot.show()

#        
    #print("Number of clusters: {0}".format(max(cluster_assignments)))
    
     
      
