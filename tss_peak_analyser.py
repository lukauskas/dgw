import pysam
import numpy
from matplotlib import pyplot

KNOWNGENES_FILE = '../data/knownGenes'

def read_tss(knowngenes_file):
    '''
    @param knowngenes_file: Filename of knownGenes
    @type knowngenes_file: str
    @rtype: set
    '''
    f = open(knowngenes_file, 'r')
    
    transcription_starting_sites = set()
    
    for line in f:
        if line[0] == '#':
            continue
        
        line = line.split('\t')
        chromosome = line[1]
        strand = line[2]
        if strand == '-':
            continue
        
        tss = int(line[3])
        
        transcription_starting_sites.add((chromosome, tss))
    
    return frozenset(transcription_starting_sites)

def process(peaks_file, alignments_file, transcription_starting_sites):
    samfile = pysam.Samfile(alignments_file, "rb")
    
    peaks_handle = open(peaks_file, "r")
    

    adjusted_n_sum = {}
    adjusted_n_count = {}
    
    adjusted_n_sum_4k = {}
    adjusted_n_count_4k = {}
   
    adjusted_n_sum_4kp = {}
    adjusted_n_count_4kp = {} 
    
    tss_counts = {}
    
    for line in peaks_handle:
        line = line.strip("\n")
        line = line.split("\t")
        
        chromosome = line[0]
        start      = int(line[1])
        end        = int(line[2])
        
        tss_locations = []
        for i in range(start, end):
            if (chromosome, i) in transcription_starting_sites:
                tss_locations.append(i)
        
        tss_count = len(tss_locations)
        
        
        try:
            tss_counts[tss_count] += 1        
        except KeyError:
            tss_counts[tss_count] = 1
            
        if tss_count != 1:
            # Skip all the peaks that have zero / more than one tss
            continue
        
        tss_loc = tss_locations[0]
        
        data = {}
        for pileupcolumn in samfile.pileup(chromosome, start, end-1):
            data[pileupcolumn.pos] = pileupcolumn.n
        
        data_arr = [ data.get(i, 0) for i in range(start, end)]
    
        for pos, n in zip(range(start, end), data_arr):
            adjusted_pos = pos - tss_loc
            
        
            try:
                adjusted_n_sum[adjusted_pos] += n
                adjusted_n_count[adjusted_pos] += 1
            except KeyError:
                adjusted_n_sum[adjusted_pos] = n
                adjusted_n_count[adjusted_pos] = 1
            
            # Get data in a similar fashion that is used in other experiments
            # that take +2000 BP -2000 BP windows of data
            if end - start <= 4001:
                try:
                    adjusted_n_sum_4k[adjusted_pos] += n
                    adjusted_n_count_4k[adjusted_pos] += 1
                except KeyError:
                    adjusted_n_sum_4k[adjusted_pos] = n
                    adjusted_n_count_4k[adjusted_pos] = 1
            else:
                try:
                    adjusted_n_sum_4kp[adjusted_pos] += n
                    adjusted_n_count_4kp[adjusted_pos] += 1
                except KeyError:
                    adjusted_n_sum_4kp[adjusted_pos] = n
                    adjusted_n_count_4kp[adjusted_pos] = 1


        #print '{0}:{1}-{2}\t{3}'.format(chromosome, start, end, len(data_arr))
        
    print '{0!r}'.format(tss_counts)
    
    x = []
    y = []
    yc = []
    
    x4k = []
    y4k = []
    y4kc = []
    
    x4kp = []
    y4kp = []
    y4kpc = []
    
    min_value = None
    max_value = None
    
    for i in sorted(adjusted_n_sum):
        x.append(i)
        average = float(adjusted_n_sum[i]) / adjusted_n_count[i]
        
        if average < min_value or min_value is None:
            min_value = average
        
        if average > max_value or max_value is None:
            max_value = average
            
        y.append(average)
        yc.append(adjusted_n_count[i])
    
    for i in sorted(adjusted_n_sum_4k):
        x4k.append(i)
        average = float(adjusted_n_sum_4k[i]) / adjusted_n_count_4k[i]
        
        if average < min_value or min_value is None:
            min_value = average
        
        if average > max_value or max_value is None:
            max_value = average
            
        y4k.append(average)
        y4kc.append(adjusted_n_count_4k[i])
    
    for i in sorted(adjusted_n_sum_4kp):
        x4kp.append(i)
        average = float(adjusted_n_sum_4kp[i]) / adjusted_n_count_4kp[i]
        
        if average < min_value or min_value is None:
            min_value = average
        
        if average > max_value or max_value is None:
            max_value = average
            
        y4kp.append(average)
        y4kpc.append(adjusted_n_count_4kp[i])
    
    pyplot.figure(1)
    pyplot.subplot(211)
    p1, = pyplot.plot(x,y)
    p2, = pyplot.plot(x4k, y4k)
    p3, = pyplot.plot(x4kp, y4kp)
    pyplot.legend([p1,p2,p3], ["All peaks", "Peaks <4k base pairs", "Peaks >4k base pairs"])
    pyplot.title("Average n per offset")
    #pyplot.vlines(0, 0, max_value)
    pyplot.subplot(212)
    p1, = pyplot.plot(x, yc)
    p2, = pyplot.plot(x4k, y4kc)
    p3, = pyplot.plot(x4kp, y4kpc)
    pyplot.title("Number of peaks containing offset")
    pyplot.legend([p1,p2], ["All peaks", "Peaks <4k base pairs", "Peaks >4k base pairs"])
    pyplot.show()
   
  
if __name__ == '__main__':
    import sys
    tss = read_tss(KNOWNGENES_FILE)
    
    process(*sys.argv[1:], transcription_starting_sites=tss)
    