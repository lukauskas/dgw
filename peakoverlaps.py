from peak import Peak
from peakanalyser import open_gzipped_file
import sys
from interval import interval

def length_of_interval(x):
    length = 0
    for pair in x:
        length += pair[1] - pair[0]
    
    return length

if __name__ == '__main__':
    filenames = sys.argv[1:]
    
    peak_intervals = {}
    
    for filename in filenames:
        #print '# Processing {0!r}'.format(filename)
        f = open_gzipped_file(filename)
        if 'broadPeak' in filename:
            line_parsing_function = Peak.parse_from_broadpeak_data_row
        elif 'narrowPeak' in filename:
            line_parsing_function = Peak.parse_from_narrowpeak_data_row
        else:
            raise Exception, 'Unknown filetype'
        
        peaks = []
        for line in f:
            peaks.append(line_parsing_function(line))
        f.close()
        
        intervals = {}
        
        intervals_store = {}
        for p in peaks:
            try:
                intervals_store[p.chromosome] |= interval([p.start, p.end])
            except KeyError:
                intervals_store[p.chromosome] = interval([p.start, p.end])
        
        peak_intervals[filename] = intervals_store

#    for (filename, peaks_store) in peak_intervals.iteritems():
#        for chromosome in sorted(peaks_store):
#            peaks = peaks_store[chromosome]
#            lenpeaks = length_of_interval(peaks)
#                
#            print '{0}\t{1}\t{2}'.format(filename, chromosome, lenpeaks)
#    
    #print '# Cross-matching'
    intersection = {}
    union        = {}
    for (filename, intervals_store) in peak_intervals.iteritems():
        if not intersection:
            for chromosome, item in intervals_store.iteritems():
                intersection[chromosome] = interval(item)
                union[chromosome]        = interval(item)
            continue
        
        for (chromosome, peaks) in intervals_store.iteritems():
            if chromosome not in intersection:
                continue
            
            intersection[chromosome] &= peaks
            union[chromosome]        |= peaks
            
    
    #print '# Calculating lengths'
    for chromosome in sorted(intersection):
        intersected_length = length_of_interval(intersection[chromosome])
        union_length = length_of_interval(union[chromosome])
        
        print "{0}\t{1}\t{2}\t{3}".format(chromosome, intersected_length, union_length, float(intersected_length) / union_length)    
        
    