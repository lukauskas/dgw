import sys
import gzip
from itertools import izip
from peak import Peak

def open_gzipped_file(filename):
    f = gzip.open(filename, 'rb')
    return f



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
        line_parsing_function = Peak.parse_from_broadpeak_data_row
    elif 'narrowPeak' in filename:
        line_parsing_function = Peak.parse_from_narrowpeak_data_row
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
