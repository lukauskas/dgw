import sys
import gzip
from itertools import izip

def open_gzipped_file(filename):
    f = gzip.open(filename, 'rb')
    return f

class Peak(object):
    chromosome = None
    _start      = None
    _end        = None

    _name       = None
    _score      = None
    strand     = None

    _signal_value = None
    _p_value      = None
    _q_value      = None

    peak         = None

    def __init__(self, chromosome, start, end, name, score, strand, signal_value, p_value, q_value, peak=None):
        self.chromosome = chromosome
        self.start = start
        self.end   = end
        self.name  = name
        self.score = score
        self.strand = strand
        self.signal_value = signal_value
        self.p_value = p_value
        self.q_value = q_value
        self.peak    = peak

    @property
    def width(self):
        return self.end - self.start

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, value):
        self._start = int(value)

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, value):
        self._end = int(value)
   
    @property
    def signal_value(self):
        return self._signal_value

    @signal_value.setter
    def signal_value(self, value):
        self._signal_value = float(value)

    @property
    def score(self):
        return self._score

    @score.setter
    def score(self, value):
        self._score = float(value)
        
    @property
    def p_value(self):
        return self._p_value

    @p_value.setter
    def p_value(self, value):
        self._p_value = float(value)
    
          
    @property
    def q_value(self):
        return self._q_value

    @q_value.setter
    def q_value(self, value):
        self._q_value = float(value)
        
    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if value == '.':
            self._name = None
        else:
            self._name = value

def broadpeak_parse(line):
    '''
    ENCODE broadPeak: Broad Peaks (or Regions) format    
     
This format is used to provide called regions of signal enrichment based on pooled, normalized (interpreted) data. It is a BED 6+3 format.

chrom - Name of the chromosome (or contig, scaffold, etc.).
chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
name - Name given to a region (preferably unique). Use '.' if no name is assigned.
score - Indicates how dark the peak will be displayed in the browser (0-1000).
strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
signalValue - Measurement of overall (usually, average) enrichment for the region.
pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
    '''
    line = line.split('\t')
    PARAM_ORDER = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 'signal_value', 'p_value', 'q_value']

    params = dict(izip(PARAM_ORDER, line))

    peak = Peak(**params)
    
    return peak

def narrowpeak_parse(line):
    line = line.split('\t')
    PARAM_ORDER = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 'signal_value', 'p_value', 'q_value', 'peak']

    params = dict(izip(PARAM_ORDER, line))

    return Peak(*line)

def bin_peaks(peaks):
    bins = {}

    for p in peaks:
        try:
            bins[p.width].append(p)
        except KeyError:
            bins[p.width] = [p]

    return bins

def mean_and_var(data):
    mu = sum(data) / len(data)
    
    varsum = 0
    for d in data:
        varsum += (d-mu) * (d - mu)
    
    var = varsum / len(data)
    
    return (mu, var)
        
if __name__ == '__main__':
    filename = sys.argv[1]
  
    f = open_gzipped_file(filename)
    if 'broadPeak' in filename:
        line_parsing_function = broadpeak_parse
    elif 'narrowPeak' in filename:
        line_parsing_function = narrowpeak_parse
    else:
        raise Exception, 'Unknown filetype'

    peaks = []
    for line in f:
        peaks.append(line_parsing_function(line))
    
    f.close()
    bins = bin_peaks(peaks)

    print 'Bin Width\tBin Size\tMean Score\tVar Score\tMean Signal\tVar Signal\tMean p-value\tVar p-value\tMean q-value\tVar q-value'
    for width in sorted(bins):
        peaks_in_the_bin = bins[width]
        
        score_mu, score_var = mean_and_var(map(lambda x: x.score, peaks_in_the_bin))        
        signal_mu, signal_var = mean_and_var(map(lambda x: x.signal_value, peaks_in_the_bin))
        p_value_mu, p_value_var = mean_and_var(map(lambda x: x.p_value, peaks_in_the_bin))
        q_value_mu, q_value_var = mean_and_var(map(lambda x: x.q_value, peaks_in_the_bin))
    

        print '{bin_width}\t{bin_size}\t{score_mu}\t{score_var}\t' \
              '{signal_mu}\t{signal_var}\t{p_value_mu}\t{p_value_var}\t' \
              '{q_value_mu}\t{q_value_var}\t' \
              .format(bin_width=width, bin_size=len(peaks_in_the_bin), 
                      score_mu=score_mu, score_var=score_var,
                      signal_mu=signal_mu, signal_var=signal_var,
                      p_value_mu=p_value_mu, p_value_var=p_value_var,
                      q_value_mu=q_value_mu, q_value_var=q_value_var)
