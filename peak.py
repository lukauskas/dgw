from itertools import izip

class Peak(object):
    table_name = None
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

    def __init__(self, table_name, chromosome, start, end, name, score, strand, signal_value, p_value, q_value, peak=None):
        self.table_name = table_name
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

    def __repr__(self):
        return '{name}({table_name}, {chromosome}, {start}, {end}, {score}, {strand}, {signal_value}, {p_value}, {q_value}, {peak})' \
            .format(name=self.__class__.__name__, track_name=self.table_name,
                  chromosome=self.chromosome, 
                  start=self.start, end=self.end, score=self.score,
                  strand=self.strand, signal_value=self.signal_value, p_value=self.p_value, q_value=self.q_value, peak=self.peak)
    
    def __add__(self, other):
        if self.chromosome != other.chromosome:
            raise ValueError, 'Cannot add peaks - chromosomes do not match'
        
        self.start = min(self.start, other.start)
        self.end   = max(self.end, other.end)
        
        # TODO: not sure what to do with other values
        return self
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
                
    @classmethod
    def parse_from_broadpeak_data_row(cls, track_name, row):
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
        
        line = row.split('\t')
        PARAM_ORDER = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 'signal_value', 'p_value', 'q_value']
    
        params = dict(izip(PARAM_ORDER, line))
    
        peak = cls(track_name, **params)
        
        return peak

    @classmethod
    def parse_from_narrowpeak_data_row(cls, track_name, row):
        line = row.split('\t')
        PARAM_ORDER = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 'signal_value', 'p_value', 'q_value', 'peak']
    
        params = dict(izip(PARAM_ORDER, line))
    
        return cls(track_name, **params)