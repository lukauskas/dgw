TSS_FILE = '../data/knownGenes'

def read_tss():
    f = open(TSS_FILE, 'r')
    
    transcription_starting_sites = set()
    
    for line in f:
        if line[0] == '#':
            continue
        
        line = line.split('\t')
        chromosome = line[1]
        tss = int(line[3])
        
        transcription_starting_sites.add((chromosome, tss))
    
    return transcription_starting_sites

def process_peak_signal(tsses, filename):
    tsses = sorted(tsses)
    peak_signals_file = open(filename, 'r')
    
    for line in peak_signals_file:
        row = line.split('\t')
        location = row[0]
        location = location.split(':')
        
        chromosome = location[0]
        
        pos = location[1].split('-')
        start = int(pos[0])
        end   = int(pos[1])
        
        tsses_for_row = []
        for tss_chrom, tss_loc in tsses:
            if tss_chrom < chromosome:
                continue
            elif tss_chrom > chromosome:
                break
            
            if tss_loc < start:
                continue
            elif tss_loc > end:
                break
            
            if tss_loc >= start and tss_loc <= end:
                tsses_for_row.append(tss_loc)
        
        if tsses_for_row:
            print '{0} | {1}'.format(line.strip('\n'), ' '.join(map(str, tsses_for_row)))
    
    peak_signals_file.close()
        
    
if __name__ == '__main__':
    import sys
    tsses = read_tss()
    
    process_peak_signal(tsses, sys.argv[1])
            
            