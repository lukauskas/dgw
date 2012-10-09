from peak import Peak
from peakanalyser import open_gzipped_file
import sys
from interval import interval

if __name__ == '__main__':
    filenames = sys.argv[1:]
    
    peak_intervals = {}
    
    for filename in filenames:
        print '# Processing {0!r}'.format(filename)
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
                intervals_store[p.chromosome].append(interval([p.start, p.end]))
            except KeyError:
                intervals_store[p.chromosome] = [interval([p.start, p.end])]
        
        for chromosome, l in intervals_store.iteritems():
            intervals_store[chromosome] = sorted(l)
        
        peak_intervals[filename] = intervals_store
    
    original_count = {}
    for (filename, peaks_store) in peak_intervals.iteritems():
        for chromosome in sorted(peaks_store):
            peaks = peaks_store[chromosome]
            lenpeaks = len(peaks)
            try:
                original_count[chromosome] += lenpeaks
            except KeyError:
                original_count[chromosome] = lenpeaks
                
            print '{0}\t{1}\t{2}'.format(filename, chromosome, lenpeaks)
    
    print '# Cross-matching'
    intersection = {}
    for (filename, intervals_store) in peak_intervals.iteritems():
        if not intersection:
            intersection = intervals_store
            continue
        
        for (chromosome, peaks) in intervals_store.iteritems():
            if chromosome not in intersection:
                continue
            
            peaks_in_intersection = intersection[chromosome]
            peaks_in_intersection = filter(lambda x: x in peaks, peaks_in_intersection)
            
            intersection[chromosome] = peaks_in_intersection
            
    for chromosome in sorted(intersection):
        intersected_length = len(intersection[chromosome])
        or_count = original_count[chromosome] - intersected_length
        
        print "{0}\t{1}\t{2}\t{3}".format(chromosome, intersected_length, or_count, float(intersected_length) / or_count)    
        
    