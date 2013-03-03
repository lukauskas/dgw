"""
A worker module for Dynamic Genome Warping.

Processes the data into an intermediate representation that could then be analysed.
Designed to be run on a multi-core machine with a lot of memory, e.g. a supercomputer.

"""
import argparse
import logging
from math import factorial
import os
import fastcluster
import numpy as np
import pandas as pd
import cPickle as pickle
from datetime import datetime
from dgw.cluster import HierarchicalClustering, compute_paths

from dgw.data.containers import Regions
from dgw.data.parsers import read_bam, HighestPileUpFilter
from dgw.dtw.parallel import parallel_pdist
from dgw.cli import StoreFilenameAction, StoreUniqueFilenameAction, Configuration


def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-r', '--regions', metavar='regions_of_interest.bed', action=StoreFilenameAction,
                        help='A BED file listing genome regions that will be processed')
    parser.add_argument('-d', '--datasets', metavar='dataset.bam', nargs='+', action=StoreUniqueFilenameAction,
                        help='One or more datasets to be analysed using DGW. Must be BAM files.')
    parser.add_argument('-pd', '--processed-dataset', metavar='processed_dataset.pd', action=StoreFilenameAction,
                        help='Dataset that has already been processed. E.g. from a previous run.')

    parser.add_argument('-poi', '--points-of-interest', metavar='poi.bed', action=StoreFilenameAction,
                        help='A BED file listing points of interest in the regions specified')

    parser.add_argument('-p', '--prefix', help='Prefix of the output files generated ', default='dgw')

    parser.add_argument('--truncate-regions', metavar='X', type=int, help='Only use first X rows of regions '
                                                                           'rather than the full dataset')

    parser.add_argument('-res', '--resolution', help='Read resolution', type=int, default=50)
    parser.add_argument('-ext', '--extend_to', help='Extend reads to specified length', type=int, default=200)

    parser.add_argument('--blank', action='store_const', const=True, default=False,
                        help='Do a blank run - just process the dataset.but do not calculate the pairwise distances')

    parser.add_argument('--metric', help='Local distance metric to be used in DTW',
                        choices=['sqeuclidean', 'euclidean', 'cosine'], default=None)

    parser.add_argument('-sb', '--slanted-band', metavar='k', help='Constrain DTW with slanted band of width k',
                        type=int) # TODO: assert > 0

    parser.add_argument('--normalise', const=True, default=False, action='store_const',
                        help='Normalise the DTW distances by dividing them by the length of longer sequence')

    parser.add_argument('-n', '--n-cpus', metavar='N', type=int,
                        help='Use up to N CPUs when calculating pairwise distances.'
                             ' Defaults to the maximum number available.')

    parser.add_argument('--output-raw-dataset', action='store_const', const=True, default=False,
                        help='Output raw dataset as well as normalised one')

    parser.add_argument('-mp', '--min-pileup', metavar='H', type=int, default=10,
                        help='Only cluster these regions that have at least one pileup column of H or more reads. ')

    parser.add_argument('-v', '--verbose', help='Turns on displaying of debug messages', action='store_const',
                        const=True, default=False)

    parser.add_argument('--output-pairwise-distances', action='store_const', const=True, default=False,
                        help='If set to true, DGW will output the pairwise distance matrix computed to a file.')

    return parser

#-- Actual execution of the program

def read_regions(regions_filename, truncate_regions):
    regions = Regions.from_bed(regions_filename)
    total_len = len(regions)
    print '> {0} regions of interest read'.format(total_len)

    used_len = total_len
    if truncate_regions:
        print '> Using only first {0} regions from {1!r}'.format(truncate_regions, regions_filename)
        used_len = truncate_regions
        regions = regions.head(truncate_regions)

    return regions, total_len, used_len

def read_datasets(configuration, regions):
    if configuration.args.datasets:
        if configuration.args.min_pileup:
            data_filters = [HighestPileUpFilter(configuration.args.min_pileup)]
        else:
            data_filters = []

        return read_bam(configuration.args.datasets, regions,
                        resolution=configuration.args.resolution,
                        extend_to=configuration.args.extend_to,
                        data_filters=data_filters,
                        output_removed_indices=True)
    else:
        processed_dataset_file = configuration.processed_dataset_filename
        logging.debug('Reading processed dataset {0!r}'.format(processed_dataset_file))

        processed_dataset = deserialise(processed_dataset_file)
        return processed_dataset, pd.Index([]), pd.Index([])  # Return empty indices for filtered regions

def serialise(object, filename):
    f = open(filename, 'wb')
    try:
        pickle.dump(object, f, protocol=pickle.HIGHEST_PROTOCOL)
    finally:
        f.close()

def deserialise(filename):
    f = open(filename, 'rb')
    try:
        obj = pickle.load(f)
        return obj
    finally:
        f.close()

def binomial_coefficent(n, k):
    return factorial(n) / (factorial(k) * factorial(n-k))

def main():
    # --- Argument parsing -----------------------
    parser = argument_parser()

    args = parser.parse_args()
    if (args.regions or args.datasets) and args.processed_dataset:
        parser.error('Must specify either both --regions and --datasets or --processed_dataset only.')
    elif not args.processed_dataset:
        if not args.regions or not args.datasets:
            parser.error('Must specify both --regions and --datasets')

    if args.metric is None:
        if args.processed_dataset:
            parser.error('Must provide a metric if using processed dataset')
        elif len(args.datasets) >= 2:
            print "> Defaulting to cosine distance as more than 2 datasets given"
            args.metric = 'cosine'
        else:
            print "> Defaulting to sqeuclidean distance as only one dataset given"
            args.metric = 'sqeuclidean'
    elif args.metric == 'cosine':
        if args.datasets and len(args.datasets) < 2:
            parser.error('Cannot use cosine distance with just one dataset. Choose sqeuclidean or euclidean instead.')

    if args.verbose:
        logging.root.setLevel(logging.DEBUG)

    configuration = Configuration(args)

    # --- pre-processing ------------------------
    if args.regions:
        print '> Reading regions from {0!r} ....'.format(args.regions)
        regions, total_regions, used_regions = read_regions(args.regions, args.truncate_regions)

        # --- saving of regions ----------
        print '> Saving regions to {0}'.format(configuration.parsed_regions_filename)
        serialise(regions, configuration.parsed_regions_filename)
    else:
        regions = None

    if args.points_of_interest:
        print '> Reading points of interest'
        poi = Regions.from_bed(args.points_of_interest)
    else:
        poi = None

    print '> Reading datasets ...'
    datasets, missing_regions, filtered_regions = read_datasets(configuration, regions)


    if args.datasets:
        if args.output_raw_dataset:
            print '> Saving raw dataset to {0}'.format(configuration.raw_dataset_filename)
            serialise(datasets, configuration.raw_dataset_filename)

        datasets = datasets.to_log_scale()

        if len(missing_regions) > 0:
            print "> {0} regions were not found in the dataset, they were saved to {1}".format(len(missing_regions),
                                                                                configuration.missing_regions_filename)
            serialise(missing_regions, configuration.missing_regions_filename)
        if len(filtered_regions) > 0:
            print "> {0} regions were filtered out from dataset due to --min-pileup constraint, they were saved to {1}".format(len(filtered_regions),
                                                                                           configuration.filtered_regions_filename)
            serialise(filtered_regions, configuration.missing_regions_filename)

        used_regions = len(datasets.items)
        if len(missing_regions) > 0 or len(filtered_regions) > 0:
            print '> {0} regions remaining'.format(used_regions)

        if poi:
            poi = poi.as_bins_of(regions, resolution=args.resolution)
            datasets.points_of_interest = poi

        # --- Saving of datasets -------------------
        print '> Saving datasets to {0}'.format(configuration.dataset_filename)
        # TODO: Pickle is likely to fuck-up here (not 64bit safe), and this is not strictly necessary!
        serialise(datasets, configuration.dataset_filename)
    else:
        print "> Not converting dataset to log scale as processed dataset already provided"

    if not args.blank:
        # --- actual work ---------------------------
        print '> Calculating pairwise distances (this might take a while) ...'
        if args.n_cpus is not None:
            print '> Using {0} processes'.format(args.n_cpus)

        start = datetime.now()
        dm = parallel_pdist(datasets, args.n_cpus, **configuration.dtw_kwargs)
        end = datetime.now()

        delta = end - start
        print '> Pairwise distances calculation took {0} s'.format(delta.total_seconds())

        if args.truncate_regions:
            multiplier = binomial_coefficent(total_regions, 2) / float(binomial_coefficent(args.truncate_regions, 2))
            print '> Expected calculation duration if truncate_regions was not used: {0} s'\
                   .format(delta.total_seconds() * multiplier)


        # --- Saving of the work --------------
        if configuration.pairwise_distances_filename:
            print '> Saving the pairwise distance matrix to {0!r}'.format(configuration.pairwise_distances_filename)
            np.save(configuration.pairwise_distances_filename, dm)

        # Linkage matrix
        print '> Computing linkage matrix'
        linkage = fastcluster.complete(dm)

        print '> Saving linkage matrix to {0!r}'.format(configuration.linkage_filename)
        np.save(configuration.linkage_filename, linkage)

        print '> Computing prototypes'
        # Hierarchical clustering object to compute the prototypes
        hc = HierarchicalClustering(datasets, regions, linkage, dtw_function=configuration.dtw_function)
        prototypes = hc.extract_prototypes()
        print '> Saving prototypes to {0!r}'.format(configuration.prototypes_filename)
        serialise(prototypes, configuration.prototypes_filename)

        print '> Computing warping paths'
        nodes = hc.tree_nodes_list
        paths = compute_paths(datasets, nodes, hc.num_obs, configuration.args.n_cpus, **configuration.dtw_kwargs)
        print '> Saving warping paths to {0!r}'.format(configuration.warping_paths_filename)
        serialise(paths, configuration.warping_paths_filename)

    else:
        print '> Skipping pairwise distances step because of --blank option set'

    print '> Saving configuration to {0!r}'.format(configuration.configuration_filename)
    f = open(configuration.configuration_filename, 'w')
    serialise(configuration, configuration.configuration_filename)

    print '> Done'


if __name__ == '__main__':
    main()