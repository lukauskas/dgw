'''
Created on 13 Nov 2012

@author: saulius
'''
import pysam
from peak import Peak
from itertools import imap
from tss_peak_analyser import read_tss, KNOWNGENES_FILE
from view.average import plot
from matplotlib import pyplot

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

def get_peaks(peaks_file, alignments_file, tss_set):
    
    samfile = pysam.Samfile(alignments_file, "rb")
    peaks_handle = open(peaks_file, "r")
    
    for line in peaks_handle:

        peak = parse_peak_from_bam_line(line)
        
        peak.find_interesting_points_from_set(tss_set)
        
        peak.data = get_peak_data(samfile, peak)
        
        yield peak

    samfile.close()
    peaks_handle.close()


if __name__ == '__main__':
    tss_set = read_tss(KNOWNGENES_FILE)
    
    import sys
    peaks = get_peaks(*sys.argv[1:], tss_set=tss_set)
    
    pyplot.figure(1)
    data = map(lambda x : x.data_relative_to_point(x.start), peaks)
    plot(data)

    pyplot.show()
      