from peak_distances import plot_hist
FILE = '../data/knownGenes'

def main():
    f = open(FILE, 'r')
    
    current_chrom = None
    last_entry = None
    distances = {}
    distances_list = []
    for line in f:
        # skip comments
        if line[0] == '#':
            continue 
        
        line = line.split('\t')
        chrom = line[1]
        start = int(line[3])
        end   = int(line[4])
        
        current_entry = (chrom,start,end)
        if last_entry is not None and chrom == current_chrom:
            distance = current_entry[1] - last_entry[2]
            
            try:
                distances[distance] += 1
            except KeyError:
                distances[distance] = 1
            distances_list.append(distance)
        
        last_entry = current_entry
        current_chrom = current_entry[0]
    
    f.close()
    
    
    for d in sorted(distances):
        count = distances[d]
        if count > 1:
            print "{0},{1}".format(d, count)
    
    plot_hist(distances_list, 'TSS distances', 'Distances between knowngenes', 10000)
    
if __name__ == '__main__':
    main()