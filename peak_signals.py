from filetools import process_peak_file
from peak import Peak
def process(peak_file, signal_file):
    
    interesting_intervals = {}
    
    prev_peak = None
    for peak in process_peak_file(peak_file):
        try: 
            interesting_intervals_for_chromosome = interesting_intervals[peak.chromosome]
        except KeyError:
            interesting_intervals[peak.chromosome] = {}
            interesting_intervals_for_chromosome = interesting_intervals[peak.chromosome]
        
        if prev_peak is not None and prev_peak.chromosome == peak.chromosome and peak.start <= prev_peak.end:
            peak.start = prev_peak.start
        
        interesting_intervals_for_chromosome[peak.start] = peak.end
        prev_peak = peak

    signal_file = open(signal_file, 'r')
    current_signal_start = None
    current_signal_end   = None
    current_signal = []
    current_signal_chromosome = None
    
    prev_chrom = None
    for row in signal_file:
        row = row.strip('\n')
        row = row.split('\t')
        
        chromosome = row[0]
        if prev_chrom != chromosome:
            interesting_intervals_for_chromosome = interesting_intervals[chromosome]
        prev_chrom = chromosome
        
        start      = int(row[1])
        end        = int(row[2])
        value      = float(row[3])
        
        if current_signal_start is not None:
            if start > current_signal_end:
                print '{0}:{1}-{2}\t{3}'.format(current_signal_chromosome, current_signal_start, current_signal_end, ' '.join(map(str, current_signal)))
                current_signal_start = None
                current_signal_end = None
                current_signal = []
                current_signal_chromosome = []
            else:
                if chromosome != current_signal_chromosome:
                    raise Exception, "Signal chromosome mismatch"
                
                current_signal.append(value)
        
        if current_signal_start is None:
            possible_end = None
            possible_start = None
            for i in range(start, end+1):
                try:
                    possible_end = interesting_intervals_for_chromosome[i]
                    possible_start = i
                    break
                except KeyError:
                    continue
                
            if possible_end is not None:
                current_signal_start = possible_start
                current_signal_end   = possible_end
                current_signal_chromosome = chromosome
    
                current_signal.append(value)
        
        
        

if __name__ == '__main__':
    import sys
    process(sys.argv[1], sys.argv[2])