from filetools import generate_description_table_from_file, process_peak_file, \
    find_description_file_from_filename
from matplotlib import pyplot
from peak import Peak
from peakanalyser import open_gzipped_file
import os
import sys
from peak_signal import GraphInterval
import pickle

OUTPUT_DIRECTORY = 'dist_histograms/'

def plot_hist(data, title, xlabel, bins=1000):
    
    #bincount = len(set(data))

    n, bins, patches = pyplot.hist(data, bins=bins)
    pyplot.xlabel(xlabel)
    pyplot.ylabel("Count")
    pyplot.title(title)
    #pyplot.savefig(os.path.join(OUTPUT_DIRECTORY, title + "_dist_hist.png"))
    
    pyplot.show()

def process(filename, description):
    last_interval = None
    current_chromosome = None
    
    distances = {}
    distances_list = []
    
    count = 0 
    for peak in process_peak_file(filename):   
        count += 1 
        if peak.chromosome != current_chromosome:
            current_chromosome = peak.chromosome
            last_interval = None
        
        if last_interval is None:
            # Not much to do for the first peak
            last_interval = GraphInterval.from_peak(peak)
            continue
        
        # If we're here we have last peak
        distance = peak.start - last_interval.end
        
        try:
            distances[distance] += 1
        except KeyError:
            distances[distance] = 1
        
        distances_list.append(distance)
            
        new_interval = GraphInterval.from_peak(peak)
#        if distance <= 0:
#            interval = last_interval + new_interval
#            print '{0}({1}+{2}): {3}'.format(
#                                    interval.genome_browser_str, 
#                                    last_interval.genome_browser_str,
#                                    new_interval.genome_browser_str,
#                                    interval.genome_browser_link)
        last_interval = new_interval
#    
    f = open(os.path.basename(filename) + '.disthist.csv', 'w')
    f.write('"Distance","Count"\n')
    for d in sorted(distances):
        if d <= 0:
            continue
        count = distances[d]
        f.write("{0},{1}\n".format(d, count))
                    
    f.close()
    
    plot_hist(distances_list, os.path.basename(filename) + '.disthist.csv', 'Distances between peaks', 10000)
    
def main():
    
    input_filename = sys.argv[1]
    
    if os.path.isdir(input_filename):
        filenames = map(lambda x : os.path.join(input_filename, x), os.listdir(input_filename))
    else:
        filenames = [input_filename]
    
    description_table = generate_description_table_from_file(find_description_file_from_filename(input_filename));
    
    for filename in filenames:
        if filename[-9:] == 'files.txt':
            continue
         
        description = description_table[os.path.basename(filename)]
        
        process(filename, description)  
if __name__ == '__main__':
    main()