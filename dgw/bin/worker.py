#!/usr/bin/env python
"""
A worker module for Dynamic Genome Warping.

Processes the data into an intermediate representation that could then be analysed.
Designed to be run on a multi-core machine with a lot of memory, e.g. a supercomputer.

"""
import argparse
import logging
from math import factorial
import random
import cPickle as pickle
from datetime import datetime
from multiprocessing import cpu_count

import fastcluster
import numpy as np
import pandas as pd

from dgw.cluster import HierarchicalClustering, compute_paths
from dgw.data.containers import Regions
from dgw.data.parsers import read_bam, HighestPileUpFilter
from dgw.data.parsers.pois import from_simple
from dgw.dtw.parallel import parallel_pdist
from dgw.cli import StoreFilenameAction, StoreUniqueFilenameAction, Configuration


def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    input_group = parser.add_argument_group('Input arguments')
    input_group.add_argument('-r', '--regions', metavar='regions_of_interest.bed', action=StoreFilenameAction,
                        help='A BED file listing genome regions that will be processed', required=False)

    input_group.add_argument('-d', '--datasets', metavar='dataset.bam', nargs='+', action=StoreUniqueFilenameAction,
                        help='One or more datasets to be analysed using DGW. Must be BAM files.')
    input_group.add_argument('-pd', '--processed-dataset', metavar='processed_dataset.pd', action=StoreFilenameAction,
                        help='Dataset that has already been processed. E.g. from a previous run. Can be used instead of -d')

    points_of_interest_group = parser.add_argument_group('Points of interest')

    points_of_interest_group.add_argument('-poi', '--points-of-interest', metavar='poi.bed', action=StoreFilenameAction,
                        help='A BED file listing points of interest (e.g. transcription start sites) '
                             'in the regions specified by --regions')
    points_of_interest_group.add_argument('--ignore-poi-non-overlaps', default=False, action='store_const', const=True,
                        help='If set to true, DGW will silently ignore points of interest that do not overlap with the regions')
    points_of_interest_group.add_argument('--ignore-no-poi-regions', default=False, action='store_const', const=True,
                        help='If set to true, DGW will silently ignore regions having no points of interest in them')

    dtw_parameters_group = parser.add_argument_group('DTW parameters')
    dtw_parameters_group.add_argument('-res', '--resolution', help='Read resolution', type=int, default=50)
    dtw_parameters_group.add_argument('-ext', '--extend_to', help='Extend reads to specified length', type=int, default=200)
    dtw_parameters_group.add_argument('--metric', help='Local distance metric to be used in DTW. Defaults to sqeuclidean '
                                                       'for one dataset, and cosine when more than one dataset is available',
                        choices=['sqeuclidean', 'euclidean', 'cosine'], default=None)
    dtw_parameters_group.add_argument('-sb', '--slanted-band', metavar='k',
                                     help='Constrain DTW with slanted band of width k', type=int) # TODO: assert > 0
    dtw_parameters_group.add_argument('-nln', '--no-length-normalisation', const=True, default=False, action='store_const',
                        help='Do not normalise the DTW distances by dividing them by the length of longer sequence.')
    dtw_parameters_group.add_argument('--no-dtw', action='store_const', const=True, default=False,
                        help='If this option is provided, DGW will not use Dynamic Time Warping for region alignments. '
                             'Use it to compare it with basic methods.')
    dtw_parameters_group.add_argument('--scale', action="store_const", const="True", default=False,
                        help='Scale the sequence uniformly to the length of the longer sequence before doing DTW')
    dtw_parameters_group.add_argument('-wp', '--warping-penalty', default=0, type=float,
                                      help='Warping penalty to incur on every non-diagonal path taken in the dtw matrix')

    dtw_parameters_group.add_argument('--use-strand-information', const=True, action="store_const", default=False,
                                  help='Will use the strand information provided in the BED file, if set to true.'
                                       ' Also sets --no-reverse.')
    dtw_parameters_group.add_argument('--no-reverse', const=True, default=False, action="store_const",
                                      help='Do not try to account for antisense promoters by reversing the regions')

    dgw_options_group = parser.add_argument_group('DGW options')
    dgw_options_group.add_argument('-p', '--prefix', help='Prefix of the output files generated ', default='dgw')

    dgw_options_group.add_argument('--random-sample', metavar='X', type=int, help='Only use a random sample of X regions '
                                                                       'rather than the full dataset')
    dgw_options_group.add_argument('--blank', action='store_const', const=True, default=False,
                        help='Do a blank run - just process the dataset.but do not calculate the pairwise distances')

    dgw_options_group.add_argument('-n', '--n-processes', metavar='N', type=int,
                        help='Use up to N process when calculating pairwise distances.'
                             ' Defaults to the maximum number available.')

    dgw_options_group.add_argument('--output-raw-dataset', action='store_const', const=True, default=False,
                        help='Output raw dataset into format readable by DGW')

    dgw_options_group.add_argument('-v', '--verbose', help='Turns on displaying of debug messages', action='store_const',
                        const=True, default=False)
    dgw_options_group.add_argument('--prototyping-method', default=None, choices=['psa', 'standard', 'standard-unweighted', 'mean'],
                                   help='Prototyping method to use. '
                                        'Defaults to standard if DGW is used and mean if it is not.')


    preprocessing_group = parser.add_argument_group('Preprocessing')
    preprocessing_group.add_argument('-mp', '--min-pileup', metavar='H', type=int, default=10,
                        help='Only cluster these regions that have at least one bin that contains H or more reads. ')
    preprocessing_group.add_argument('--normalise-pileups', action='store_const', const=True, default=False,
                                     help='Normalise the number of elements in the bins by dividing them by the number of reads '
                                          'falling into the bin with the largest reads.')


    preprocessing_group.add_argument('--output-pairwise-distances', action='store_const', const=True, default=False,
                        help='If provided, DGW will output the intermediate pairwise distance matrix computed.')




    preprocessing_group.add_argument('--min-bins', default=4, help='Specifies the minimal length of region in bins in order for it '
                                                      'to be processed.', type=int)
    preprocessing_group.add_argument('--max-bins', default=1000, help='Specifies the maximum length of region in bins '
                                                          'in order for it to be processed',
                        type=int)





    return parser

#-- Actual execution of the program

def read_regions(regions_filename, random_sample, resolution):
    regions = Regions.from_bed(regions_filename)
    total_len = len(regions)
    print '> {0} regions of interest read'.format(total_len)

    regions = regions.clip_to_resolution(resolution)

    used_len = total_len
    if random_sample:
        print '> Using only a random sample of {0} regions from {1!r}'.format(random_sample, regions_filename)
        used_len = random_sample
        regions = regions.ix[random.sample(regions.index, random_sample)]

    return regions, total_len, used_len

def read_datasets(args, regions):
    if args.datasets:
        data_filters = []
        if args.min_pileup:
            data_filters.append(HighestPileUpFilter(args.min_pileup))

        return read_bam(args.datasets, regions,
                        resolution=args.resolution,
                        extend_to=args.extend_to,
                        data_filters=data_filters,
                        output_removed_indices=True,
                        reverse_negative_strand_regions=args.use_strand_information)
    else:
        processed_dataset_file = args.processed_dataset
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
    if args.datasets and args.processed_dataset:
        parser.error('Must specify either --dataset or --processed_dataset only.')
    elif not args.processed_dataset:
        if not args.regions or not args.datasets:
            parser.error('Must specify both --regions and --dataset')

    if args.metric is None:
        if args.processed_dataset:
            parser.error('Must provide a metric if using processed dataset')
        elif len(args.datasets) >= 2:
            print "> Defaulting to cosine distance as more than 2 dataset given"
            args.metric = 'cosine'
        else:
            print "> Defaulting to sqeuclidean distance as only one dataset given"
            args.metric = 'sqeuclidean'
    elif args.metric == 'cosine':
        if args.datasets and len(args.datasets) < 2:
            parser.error('Cannot use cosine distance with just one dataset. Choose sqeuclidean or euclidean instead.')

    if args.no_dtw:
        # That's what no-dtw actually does
        args.slanted_band = 0
        args.scale = True
        if args.prototyping_method is None:
            args.prototyping_method = 'mean'
    else:
        if args.prototyping_method is None:
            args.prototyping_method = 'standard'

    if args.verbose:
        logging.root.setLevel(logging.DEBUG)

    # Disable trying to reverse regions if strand information given
    if args.use_strand_information:
        args.no_reverse = True

    configuration = Configuration(args)

    # --- pre-processing ------------------------
    if args.regions:
        print '> Reading regions from {0!r} ....'.format(args.regions)
        regions, total_regions, used_regions = read_regions(args.regions, args.random_sample, args.resolution)
        if args.use_strand_information and not regions.has_strand_data():
            logging.debug('Parsed columns: {0}'.format(regions.columns))
            parser.error('--use-strand-information is set but the input BED file has no strand information.')

        too_short_regions = (regions.lengths / args.resolution) < args.min_bins  # Set the threshold to 4 bins
        too_short_regions = regions.ix[too_short_regions[too_short_regions].index]
        if len(too_short_regions) > 0:
            print '> {0} regions have their length shorter than {1} bins. Saving them to {2!r} as they won\'t be processed'\
                .format(len(too_short_regions), args.min_bins, configuration.too_short_regions_filename)
            too_short_regions.to_bed(configuration.too_short_regions_filename)

            regions = regions.ix[regions.index - too_short_regions.index]

        if args.max_bins:
            too_long_regions = (regions.lengths / args.resolution) >= args.max_bins
            too_long_regions = regions.ix[too_long_regions[too_long_regions].index]

            if len(too_long_regions) > 0:
                print '> {0} regions have their length longer than {1} bins. ' \
                      'Saving them to {2!r} as they won\'t be processed due to --max-bins constraint'\
                      .format(len(too_long_regions), args.max_bins, configuration.too_long_regions_filename)
                too_long_regions.to_bed(configuration.too_long_regions_filename)

                regions = regions.ix[regions.index - too_long_regions.index]

        print '> {0} regions remain'.format(len(regions))
    else:
        regions = None

    if args.points_of_interest:
        print '> Reading points of interest'

        poi_file = args.points_of_interest
        try:
            poi = from_simple(poi_file, regions, resolution=configuration.resolution)
        except ValueError:
            poi = Regions.from_bed(poi_file)
            poi = poi.as_bins_of(regions, resolution=configuration.resolution,
                                 ignore_non_overlaps=args.ignore_poi_non_overlaps,
                                 account_for_strand_information=configuration.use_strand_information)

        if not poi:
            raise Exception(
                'POI file provided, but no POIs were parsed from {}. Try using dgw-overlaps2poi'.format(
                    poi_file))
    else:
        poi = None

    print '> Reading dataset ...'
    dataset, missing_regions, filtered_regions = read_datasets(args, regions)

    if args.datasets:

        if poi:
            dataset.add_points_of_interest(poi, name=args.points_of_interest)

            if args.ignore_no_poi_regions:
                poi_dataset = dataset.drop_no_pois()

                if len(poi_dataset) != len(dataset):
                    dropped_regions = regions.ix[dataset.items - poi_dataset.items]
                    print '> {0} regions were removed as they have no POI data with them ' \
                          'and --ignore-no-poi-regions was set'.format(len(dropped_regions))
                    print '> Saving them to {0!r}'.format(configuration.no_poi_regions_filename)
                    dropped_regions.to_bed(configuration.no_poi_regions_filename)
                    dataset = poi_dataset
                    del dropped_regions
                del poi_dataset

        if len(missing_regions) > 0:
            print "> {0} regions were not found in the dataset, they were saved to {1}".format(len(missing_regions),
                                                                                configuration.missing_regions_filename)
            regions.ix[missing_regions].to_bed(configuration.missing_regions_filename, track_title='DGWMissingRegions',
                                   track_description='Regions that are in input, but missing from the dataset')

        if len(filtered_regions) > 0:
            print "> {0} regions were filtered out from dataset due to --min-pileup constraint, they were saved to {1}".format(len(filtered_regions),
                                                                                           configuration.filtered_regions_filename)
            regions.ix[filtered_regions].to_bed(configuration.filtered_regions_filename, track_title='DGWFilteredRegions',
                                        track_description='Regions that were filtered out from the dataset')

        # Get remaining regions
        regions = regions.ix[dataset.items]
        if len(missing_regions) > 0 or len(filtered_regions) > 0:
            print '> {0} regions remaining and will be processed'.format(len(regions))


        if args.output_raw_dataset:
            print '> Serialising raw dataset to {0}'.format(configuration.raw_dataset_filename)
            serialise(dataset, configuration.raw_dataset_filename)

        dataset = dataset.to_log_scale()

        if args.normalise_pileups:
            print '> Dividing the number of reads in each bin by the maximum number of reads per region as --normalise-pileups is set'
            dataset = dataset.normalise_bin_heights()


        missing_regions = regions.index - dataset.items

        if len(missing_regions) > 0:
            print "> {0} regions were not found in the dataset, they were saved to {1}".format(len(missing_regions),
                                                                                               configuration.missing_regions_filename)
            regions.ix[missing_regions].to_bed(configuration.missing_regions_filename, track_title='DGWMissingRegions',
                                               track_description='Regions that are in input, but missing from the dataset')
    else:
        print "> Not converting dataset to log scale as processed dataset already provided"

    # --- Serialise the regions as they will be needed in explorer ----------
    if regions is not None:
        print '> Serialising regions to {0}'.format(configuration.parsed_regions_filename)
        serialise(regions, configuration.parsed_regions_filename)

    # --- Saving of dataset -------------------
    print '> Saving dataset to {0}'.format(configuration.dataset_filename)
    serialise(dataset, configuration.dataset_filename)

    if not args.blank:
        # --- actual work ---------------------------
        print '> Calculating pairwise distances (this might take a while) ...'
        if args.n_processes is not None:
            print '> Using {0} processes'.format(args.n_processes)
        else:
            args.n_processes = cpu_count()
            print '> Using all available cpu cores ({0})'.format(args.n_processes)

        if args.no_dtw:
            print '> Not using DTW as --no-dtw option is set'

        logging.debug('Running DTW with the following kwargs: {0!r}'.format(configuration.dtw_kwargs))
        start = datetime.now()
        dm = parallel_pdist(dataset, args.n_processes, **configuration.dtw_kwargs)
        end = datetime.now()

        delta = end - start
        print '> Pairwise distances calculation took {0} s'.format(delta.total_seconds())

        if args.random_sample:
            multiplier = binomial_coefficent(total_regions, 2) / float(binomial_coefficent(args.random_sample, 2))
            print '> Expected calculation duration if random-sample was not used: {0} s'\
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
        hc = HierarchicalClustering(dataset, regions, linkage, dtw_function=configuration.dtw_function,
                                    prototyping_method=configuration.prototyping_method)
        prototypes = hc.extract_prototypes()
        print '> Saving prototypes to {0!r}'.format(configuration.prototypes_filename)
        serialise(prototypes, configuration.prototypes_filename)

        print '> Computing warping paths'
        nodes = hc.tree_nodes_list
        paths = compute_paths(dataset, nodes, hc.num_obs, n_processes=args.n_processes,
                              **configuration.dtw_kwargs)
        print '> Saving warping paths to {0!r}'.format(configuration.warping_paths_filename)
        serialise(paths, configuration.warping_paths_filename)
    else:
        print '> Skipping pairwise distances step because of --blank option set'

    print '> Saving configuration to {0!r}'.format(configuration.configuration_filename)
    f = open(configuration.configuration_filename, 'w')
    try:
        configuration.to_json(f)
    finally:
        f.close()

    print '> Done'


if __name__ == '__main__':
    main()