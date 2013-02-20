"""
A worker module for Dynamic Genome Warping.

Processes the data into an intermediate representation that could then be analysed.
Designed to be run on a multi-core machine with a lot of memory, e.g. a supercomputer.

"""
import argparse
from math import factorial
import os
import numpy as np
import pandas as pd
import cPickle as pickle
from datetime import datetime

from dgw.data.containers import Regions
from dgw.data.parsers import read_bam
from dgw.dtw.parallel import parallel_pdist
from dgw.cli import StoreFilenameAction, StoreUniqueFilenameAction, Configuration

def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--regions', metavar='regions_of_interest.bed', action=StoreFilenameAction,
                        help='A BED file listing genome regions that will be processed', required=True)
    parser.add_argument('-d', '--datasets', metavar='dataset.bam', nargs='+', action=StoreUniqueFilenameAction,
                        help='One or more datasets to be analysed using DGW. Must be BAM files.', required=True)

    parser.add_argument('-p', '--prefix', help='Prefix of the output files generated '
                                               '(defaults to \'dgw\')')

    parser.add_argument('--truncate-regions', metavar='X', type=int, help='Only use first X rows of regions '
                                                                           'rather than the full dataset')

    parser.add_argument('-res', '--resolution', help='Read resolution', type=int, default=50)
    parser.add_argument('-ext', '--extend_to', help='Extend reads to specified length', type=int, default=200)

    parser.add_argument('--blank', action='store_const', const=True, default=False,
                        help='Do a blank run - just process the dataset.but do not calculate the pairwise distances')

    parser.add_argument('--metric', help='Local distance metric to be used in DTW',
                        choices=['sqeuclidean', 'euclidean', 'cosine'], default='sqeuclidean')

    parser.add_argument('-sb', '--slanted-band', metavar='k', help='Constrain DTW with slanted band of width k',
                        type=int) # TODO: assert > 0

    parser.add_argument('-n', '--n-cpus', metavar='N', type=int,
                        help='Use up to N CPUs when calculating pairwise distances.'
                             ' Defaults to the maximum number available.')

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

def read_datasets(dataset_filenames, regions, resolution, extend_to):
    return read_bam(dataset_filenames, regions, resolution, extend_to)

def serialise(object, filename):
    f = open(filename, 'w')
    pickle.dump(object, f, protocol=pickle.HIGHEST_PROTOCOL)
    f.close()

def compare_datasets_and_regions(regions, datasets):
    missing_regions = regions.regions_not_in_dataset(datasets)
    print '> {0} regions were not found in the dataset'.format(len(missing_regions))

    return missing_regions

def binomial_coefficent(n, k):
    return factorial(n) / (factorial(k) * factorial(n-k))

def main():
    # --- Argument parsing -----------------------
    parser = argument_parser()

    args = parser.parse_args()
    if args.prefix is None:
        args.prefix = 'dgw'

    configuration = Configuration(args)

    # --- pre-processing ------------------------
    print '> Reading regions from {0!r} ....'.format(args.regions)
    regions, total_regions, used_regions = read_regions(args.regions, args.truncate_regions)

    # --- saving of regions ----------
    print '> Saving regions to {0}'.format(configuration.parsed_regions_filename)
    serialise(regions, configuration.parsed_regions_filename)

    print '> Reading datasets ...'
    datasets = read_datasets(args.datasets, regions, args.resolution, args.extend_to)
    datasets = datasets.to_log_scale()

    missing_regions = compare_datasets_and_regions(regions, datasets)
    serialise(missing_regions, configuration.missing_regions_filename)

    if len(missing_regions) > 0:
        print "> There were {0} regions that are in input file, but not in dataset.".format(len(missing_regions))
        print "> These regions were saved to are saved to file {0}".format(configuration.missing_regions_filename)

        used_regions -= missing_regions

    # --- Saving of datasets -------------------
    print '> Saving datasets to {0}'.format(configuration.dataset_filename)
    # TODO: Pickle is likely to fuck-up here (not 64bit safe), and this is not strictly necessary!
    serialise(datasets, configuration.dataset_filename)

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
        print '> Saving the pairwise distance matrix to {0!r}'.format(configuration.pairwise_distances_filename)
        np.save(configuration.pairwise_distances_filename, dm)
    else:
        print '> Skipping pairwise distances step because of --blank option set'

    print '> Saving configuration to {0!r}'.format(configuration.configuration_filename)
    f = open(configuration.configuration_filename, 'w')
    serialise(configuration, configuration.configuration_filename)

    print '> Done'


if __name__ == '__main__':
    main()