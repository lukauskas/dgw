import matplotlib.pyplot as plt
from numpy.core.numeric import arange

GENOME_BROWSER_EMPTY_SESSION_URI = 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=sauliusl&hgS_otherUserSessionName=hg19_empty'

class Signal(object):
    
    chromosome = None
    data = None
    
    def __init__(self, chromosome):
        self.chromosome = chromosome
        self.data = []
        
    def add_data_point(self, start, end, signal_value):
        start = int(start)
        end   = int(end)
        signal_value = float(signal_value)

        self.data.append((start, end, signal_value))
        
def plot_signal(signal):
    data = signal.data
    points  = [ x[0]  for x in data ] 
    heights = [ x[2] for x in data ]
    plt.plot(points, heights)
    plt.yticks(arange(40))
    plt.show()

class GraphInterval(object):
    
    track_name  = None
    chromosome = None
    start = None
    end   = None
    
    def __init__(self, track_name, chromosome, start, end):
        self.track_name = track_name
        self.chromosome = chromosome
        self.start      = int(start)
        self.end        = int(end)
    
    def contains(self, chromosome, start, end):
        return self.chromosome == chromosome and self.start <= int(start) and self.end > int(start)
    
    @property
    def genome_browser_str(self):
        return '{0}:{1}-{2}'.format(self.chromosome, self.start, self.end)
    
    @property
    def genome_browser_link(self):
        peak_track = self.track_name + 'Pk'
        signal_track = self.track_name + 'Sig'
        
        return '{base_uri}&position={position}&{peak_track}=pack&{signal_track}=full' \
                    .format(base_uri=GENOME_BROWSER_EMPTY_SESSION_URI,
                            position=self.genome_browser_str,
                            peak_track=peak_track, signal_track=signal_track)
        
        
    def __repr__(self):
        return '<GraphInterval {0}>'.format(self.genome_browser_str)
    
    @classmethod
    def from_peak(cls, peak):
        track_name = peak.table_name[:-2]
        return cls(track_name, peak.chromosome, peak.start, peak.end)
    
    def __add__(self, other):
        if self.track_name != other.track_name:
            raise TypeError, 'Track names do not match'
        if self.chromosome != other.chromosome:
            raise TypeError, 'Chromosomes do not match'
        
        return GraphInterval(self.track_name, self.chromosome, min(self.start, other.start), max(self.end, other.end))
    
    
   
def find_signals_for_intervals(bed_graph_file, intervals):
    
    f = open(bed_graph_file, 'r')
    intervals = sorted(intervals)
    
    signals = []
    chromosomes_of_interest = set()
    
    for interval in intervals:
        chromosomes_of_interest.add(interval.chromosome)
        signals.append(Signal(interval.chromosome))
          
    for line in f:
        line = line.split('\t')
     
        # Weed out wrong chromosomes
        line_chromosome = line[0]
        if line_chromosome not in chromosomes_of_interest:
            continue
        
        line_start = int(line[1])
        line_end   = int(line[2])
        line_value = float(line[3])
        
        for (interval, signal) in zip (intervals, signals):
            if interval.contains(line_chromosome, line_start, line_end):
                signal.add_data_point(line_start, line_end, line_value)
     
    f.close()
    
    return zip(intervals, signals)