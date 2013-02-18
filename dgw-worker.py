"""
A worker module for Dynamic Genome Warping.

Processes the data into an intermediate representation that could then be analysed.
Designed to be run on a multi-core machine with a lot of memory, e.g. a supercomputer.

"""
import argparse
import os
import numpy as np
import pandas as pd
import cPickle as pickle

from dgw.data.containers import Regions
from dgw.data.parsers import read_bam
from dgw.dtw.parallel import parallel_pdist

class StoreUniqueFilenameAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        seen = set()
        for value in values:
            if value in seen:
                raise argparse.ArgumentError(self, 'Dataset {0!r} is twice in the list'.format(value))
            if not os.path.isfile(value):
                raise argparse.ArgumentError(self, 'File not found {0!r}'.format(value))

            seen.add(value)

        setattr(namespace, self.dest, values)

class StoreFilenameAction(argparse.Action):

    def __call__(self, parser, namespace, value, option_string=None):
        if not os.path.isfile(value):
            raise argparse.ArgumentError(self, 'File not found {0!r}'.format(value))

        setattr(namespace, self.dest, value)


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
                        choices=['sqeuclidean', 'euclidean'], default='sqeuclidean')

    return parser

#-- Actual execution of the program

def read_regions(regions_filename, truncate_regions):
    regions = Regions.from_bed(regions_filename)

    if truncate_regions:
        print 'Using only first {0} regions from {1!r}'.format(truncate_regions, regions_filename)
        regions = regions.head(truncate_regions)

    return regions

def read_datasets(dataset_filenames, regions, resolution, extend_to):
    return read_bam(dataset_filenames, regions, resolution, extend_to)

class Configuration(object):
    _args = None

    def __init__(self, args):
        self._args = args

    @property
    def args(self):
        return self._args

    @property
    def pairwise_distances_filename(self):
        return '{0}_pairwise_distances.npy'.format(self.args.prefix)
    @property
    def configuration_filename(self):
        return '{0}_config.pickle'.format(self.args.prefix)

    @property
    def parsed_regions_filename(self):
        return '{0}_regions.pd'.format(self.args.prefix)

    @property
    def dataset_filename(self):
        return '{0}_datasets.pd'.format(self.args.prefix)

def serialise(object, filename):
    f = open(filename, 'w')
    pickle.dump(object, f, protocol=pickle.HIGHEST_PROTOCOL)
    f.close()

def main():
    # --- Argument parsing -----------------------
    parser = argument_parser()
    args = parser.parse_args()
    if args.prefix is None:
        args.prefix = 'dgw'

    configuration = Configuration(args)

    # --- pre-processing ------------------------
    print '> Reading regions from {0!r} ....'.format(args.regions)
    regions = read_regions(args.regions, args.truncate_regions)

    # --- saving of regions ----------
    print '> Saving regions to {0}'.format(configuration.parsed_regions_filename)
    serialise(regions, configuration.parsed_regions_filename)

    print '> Reading datasets ...'
    datasets = read_datasets(args.datasets, regions, args.resolution, args.extend_to)
    datasets = datasets.to_log_scale()

    # --- Saving of datasets -------------------
    print '> Saving datasets to {0}'.format(configuration.dataset_filename)
    # TODO: Pickle is likely to fuck-up here (not 64bit safe), and this is not strictly necessary!
    serialise(datasets, configuration.dataset_filename)

    if not args.blank:
        # --- actual work ---------------------------
        print '> Calculating pairwise distances (this might take a while) ...'
        dm = parallel_pdist(datasets, metric=args.metric)

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