import pickle
import sys

from peak_signal import *

def main(filename, signal_filename):
    f = open(filename, 'rb')
    intervals = pickle.load(f)
    f.close()
    
    overlaps = intervals['overlap']
    
    signals = find_signals_for_intervals(signal_filename, overlaps)
    
    for interval, signal in signals:
        print interval
        plot_signal(signal)
        
    
if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])