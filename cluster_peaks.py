'''
Created on 13 Nov 2012

@author: saulius
'''
import pysam
from peak import Peak
from itertools import izip
from tss_peak_analyser import read_tss, KNOWNGENES_FILE
from view.average import plot
from matplotlib import pyplot
from cluster.distance import dtw_distance_matrix
from cluster.cluster import hierarchical_cluster
import math
from cache import cached
import os

def parse_peak_from_bam_line(line):
    '''
    @rtype: Peak
    '''
    line = line.strip("\n")
    line = line.split("\t")
    
    chromosome = line[0]
    start      = line[1]
    end        = line[2]
    
    return Peak(chromosome, start, end)

def get_peak_data(samfile, peak):
    '''
    @param samfile: samfile to query
    @type samfile: pysam.Samfile
    @param peak: peak to get the data for
    @type peak: Peak
    @rtype: list
    '''
    
    data = {}
    
    # Note this has to be peak.end - 1 
    # because samfile.pileup treats ends inclusively
    for pileupcolumn in samfile.pileup(peak.chromosome, 
                                       peak.start, 
                                       peak.end-1):
        
        data[pileupcolumn.pos] = int(pileupcolumn.n)
    
    data_arr = [ data.get(i, 0) 
                 for i in range(peak.start, peak.end)]

    return data_arr    

@cached
def get_peaks(peaks_file, alignments_file, tss_set):
    
    samfile = pysam.Samfile(alignments_file, "rb")
    peaks_handle = open(peaks_file, "r")
    
    peaks = []
    for line in peaks_handle:

        peak = parse_peak_from_bam_line(line)
        
        peak.find_interesting_points_from_set(tss_set)
        
        peak.data = get_peak_data(samfile, peak)
        
        peaks.append(peak)
        

    samfile.close()
    peaks_handle.close()
    
    return peaks

if __name__ == '__main__':
    tss_set = read_tss(KNOWNGENES_FILE)
    
    import sys
    
    print("Obtaining data")
    peaks = list(get_peaks(*sys.argv[1:], tss_set=tss_set))
    
    data = map(lambda x : x.data, peaks)
    N = len(data)
    print("Calculating distance")
    
    THRESHOLD = 500
    distance_func = cached(dtw_distance_matrix, os.path.basename(sys.argv[1]))
    cluster_assignments = hierarchical_cluster(data, distance_func, THRESHOLD)
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
        cluster_data = [ p.data_relative_to_start() for c,p in izip(cluster_assignments, peaks) if c == i]
        print("Peak {0}: {1}".format(i, len(cluster_data)))
        plot(cluster_data)
    
    pyplot.show()
        
    print("Number of clusters: {0}".format(max(cluster_assignments)))
    
     
      