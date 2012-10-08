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
    score      = None
    strand     = None

    _signal_value = None
    p_value      = None
    q_value      = None

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
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if value == '.':
            self._name = None
        else:
            self._name = value

def broadpeak_parse(line):
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

def mean(data):
    return sum(data) / len(data)

def variance(data):
    mu = mean(data)
    
    varsum = 0
    for d in data:
        varsum += (d-mu) * (d - mu)

    return varsum / len(data)
        
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

    for width in sorted(bins):
        peaks_in_the_bin = bins[width]
        peak_signals = map(lambda x: x.signal_value, peaks_in_the_bin)
        mu = mean(peak_signals)
        var = variance(peak_signals)

        print '{0}\t{1}\t{2}\t{3}'.format(width, len(peaks_in_the_bin), mu, var)
