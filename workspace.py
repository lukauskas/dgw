from dgw.data.parsers import genes

__author__ = 'saulius'

KNOWN_GENES = '../data/knownGenes'
K562_H3K4ME3_REP1 = '../data/interesting/broad/K562/wgEncodeBroadHistoneK562H3k4me3StdAlnRep1.bam'
K562_H3K4ME3_REP2 = '../data/interesting/broad/K562/wgEncodeBroadHistoneK562H3k4me3StdAlnRep2.bam'
K562_CONTROL_REP1 = '../data/interesting/broad/K562/wgEncodeBroadHistoneK562ControlStdAlnRep1.bam'
MACS_MACS_H3K4ME3_REP1 = '../data/interesting/broad/K562/K56H3k4me3Rep1_peaks.bed'
CPG_ISLANDS = '../data/cpg_islands'
import dgw.data.parsers.cpgs
RESOLUTION = 50

import pandas as pd
import matplotlib.pyplot as plt
import dgw.data.visualisation.heatmap as heatmap
import webbrowser
import numpy as np
import dgw.dtw.transformations

CLEAN_VALID_GENE_REGIONS_FILENAME = 'clean_valid_gene_regions.pandas'

def overlapping_regions(region, data, window=2000):
    chr_regions = data[data.chromosome==region['chromosome']]
    overlapping_regions = chr_regions[chr_regions.start.between(region['start']-window, region['end'])]
    return overlapping_regions.ix[overlapping_regions.index - [region.name]]

def overlapping_regions_count(regions, data):

    series = []
    ixs    = []

    for ix, region in regions.iterrows():
        overlaps = overlapping_regions(region, data)
        series.append(len(overlaps))
        ixs.append(ix)

    return pd.Series(series, index=ixs)

def xticks_rel_to_tss(resolution, window):
    (xticks, xtick_labels) = plt.xticks()
    xtick_labels = [ int(x*resolution - window) for x in xticks ]
    plt.gca().set_xticklabels(xtick_labels)
    plt.gca().set_xlabel('Position relative to TSS')

#def plot_cluster(sub_data, cluster_id, plots=['heatmap', 'mean'], resolution=None, window=None):
#    plt.figure()
#
#    plot_count = len(plots)
#    if plot_count == 0:
#        raise ValueError('Nothing to plot - plots is empty')
#
#    for j, plot in enumerate(plots):
#        plt.subplot(1, plot_count, j+1)
#        if plot == 'heatmap':
#            heatmap.plot(sub_data)
#            if resolution is not None and window is not None:
#                xticks_rel_to_tss(resolution, window)
#        elif plot == 'mean':
#            sub_data.mean().plot()
#            if resolution is not None and window is not None:
#                xticks_rel_to_tss(resolution, window)
#        elif plot == 'dba':
#            dba, _ = dtw.dba(sub_data.values, sub_data.values[0])
#            plt.plot(dba)
#            if resolution is not None and window is not None:
#                xticks_rel_to_tss(resolution, window)
#
#    plt.suptitle('Cluster {0} ({1} items)'.format(cluster_id, len(sub_data)))
#
#
#
#def plot_clusters(clusters, data, plots=['heatmap', 'mean'], resolution=None, window=None, filter=None):
#
#    VALID_PLOTS = set(['heatmap', 'mean', 'dba'])
#    for plot in plots:
#        if plot not in VALID_PLOTS:
#            raise ValueError('Unknown plot {0}'.format(plot))
#
#    nclusters = clusters.max()
#
#    if (resolution is not None and window is None) or (resolution is None and window is not None):
#        raise ValueError('Both resolution and window should be provided, not just one of them')
#
#    for i in range(1, nclusters+1):
#
#        sub_data = data[clusters==i]
#
#        if filter is not None:
#            if not filter(sub_data):
#                continue
#
#        plot_cluster(sub_data, i, plots=plots, resolution=resolution, window=window)


def plot_cluster(cluster):
    plt.figure()
    plt.subplot(2,2,1)
    heatmap.plot(cluster.items)
    plt.title('Items in cluster')
    plt.subplot(2,2,2)
    plt.plot(cluster.prototype())
    plt.title('Prototype')
    plt.subplot(2,2,3)
    base = cluster.prototype()
    projections = dgw.dtw.transformations.dtw_projection_multi(cluster.items, base)
    heatmap.plot(projections)
    plt.title('Items projected onto prototype')
    plt.subplot(2,2,4)
    projections.mean().plot()
    plt.title('Mean projections on prototype')


def genome_browser_link(region):
    return 'http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=sauliusl&hgS_otherUserSessionName=hg19_k562-h3k4me3&position={0}:{1}-{2}'.format(region['chromosome'], region['start'], region['end'])
def open_in_genome_browser(region):
    webbrowser.open(genome_browser_link(region))

if __name__ == '__main__':
    print '> Initialising..'

    known_genes = genes.read_encode_known_genes(KNOWN_GENES)

    import sys
    if sys.argv[1] == '--tss':
        # Initialise
        clean_valid_gene_regions = pd.load(CLEAN_VALID_GENE_REGIONS_FILENAME)
        peak_data = pd.load('peak_data_50-200.pandas')
#        peak_data = read_bam(K562_H3K4ME3_REP1, clean_valid_gene_regions, resolution=50)
        log_peak_data = (peak_data + 1).apply(np.log)
        known_genes = known_genes.ix[peak_data.index]
    elif sys.argv[1] == '--macs':
        regions =  read_bed(MACS_MACS_H3K4ME3_REP1, resolution=RESOLUTION)
        peak_data = pd.load('macs_peak_data_50.pandas')
        #peak_data = read_bam(K562_H3K4ME3_REP1, regions, resolution=RESOLUTION, extend_to=200)
        #norm_peak_data = peak_data.div(peak_data.T.sum(), axis='index')
    elif sys.argv[1] == '--cpgs':
        cpgs = dgw.data.parsers.cpgs.read_cpgs(CPG_ISLANDS)
        cpg_regions = pd.DataFrame({'chromosome' : cpgs['chromosome'], 'start' : cpgs['start'] - 2000, 'end' : cpgs['end'] - 2000})
#        cpg_data = pd.load('cpg_data.pandas')
#        cpg_data_log = (cpg_data + 1).apply(np.log)
#        dm = np.load('dm_cpg_data_log_50.npy')



