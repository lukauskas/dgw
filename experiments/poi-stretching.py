import argparse
from collections import defaultdict
import os
from dgw import HighestPileUpFilter, read_bam, Regions
from dgw.cli import StoreFilenameAction, StoreUniqueFilenameAction
import cPickle as pickle
import random
import matplotlib.pyplot as plt
from dgw.data.containers import AlignmentsData
from dgw.data.parsers.pois import from_simple
from dgw.dtw import uniform_scaling_to_length, reverse_sequence
from dgw.evaluation.resampling import mutate_sequence
import numpy as np
import pandas as pd

def load_from_pickle(filename):
    f = open(filename)
    try:
        return pickle.load(f)
    finally:
        f.close()

def read_datasets(args, regions):
    data_filters = []
    if args.min_pileup:
        data_filters.append(HighestPileUpFilter(args.min_pileup))

    return read_bam(args.datasets, regions,
                    resolution=args.resolution,
                    extend_to=args.extend_to,
                    data_filters=data_filters,
                    output_removed_indices=False,
                    reverse_negative_strand_regions=args.use_strand_information)

def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--datasets', metavar='dataset.bam', nargs='+', action=StoreUniqueFilenameAction,
                             help='One or more datasets to be used in the rewarped dataset. Must be BAM files.',
                             required=True)
    parser.add_argument('-r', '--regions', metavar='regions.bed', action=StoreFilenameAction,
                        help='Regions of interest to use for rewarping', required=True)
    parser.add_argument('-poi', '--points-of-interest', metavar='poi.poi', action=StoreFilenameAction,
                        help='Points of interest that will be stretched. Only continuous BED points of interest allowed',
                        required=True)
    parser.add_argument('-O', '--output', metavar='processed_dataset.pd', required=True,
                        help='Output file to store the processed dataset at')
    parser.add_argument('-OR', '--output-regions', metavar='regions.bed', required=True,
                        help='File in which the original regions that were used to generate the dataset will be output')
    parser.add_argument('-n', '--number-of-regions', type=int, help='Specified how many regions to rewarp to generate the'
                                                                    ' new dataset.', required=True)
    parser.add_argument('-s', '--size', type=int, help='Size of the newly created dataset', required=True)

    parser.add_argument('-mp', '--min-pileup', metavar='H', type=int, default=10,
                                     help='Only cluster these regions that have at least one bin that contains H or more reads. ')

    parser.add_argument('-res', '--resolution', help='Read resolution', type=int, default=50)
    parser.add_argument('-ext', '--extend_to', help='Extend reads to specified length', type=int, default=200)

    parser.add_argument('--max-stretch-length', metavar='N', type=int, help='Max length to stretch the POI regions to (in bins)',
                        default=20)
    parser.add_argument('--no-reverse', action='store_const', default=False, const=True,
                        help='Do not randomly reverse the regions')
    parser.add_argument('--mutation-probability', type=float, default=0.01, help='Probability of mutation of the value '
                                                                                 'of some bin in the region')
    parser.add_argument('--mutation-scale', type=float, default=0.2, help='Scale of the mutations. Value of 0.2 means'
                                                                                'that the items could mutate from -20% to 20%')
    parser.add_argument('--ignore-poi-non-overlaps', default=False, action='store_const', const=True,
                         help='If set to true, silently ignore points of interest that do not overlap with the regions')

    parser.add_argument('--use-strand-information', const=True, action="store_const",
                                      help='Will use the strand information provided in the BED file, if set to true.')

    return parser

def get_random_sample(dataset, n):
    """
    Returns a random sample of the dataset of size n

    :param dataset:
    :param n:
    :return:
    """
    return dataset.ix[random.sample(dataset.items, n)]

def plot_heatmap(dataset):
    """
    Plots heatmap of dataset with POIs superimposed on it
    :param dataset:
    :return:
    """
    dataset.plot_heatmap(highlighted_points=dataset.points_of_interest, sort_by=None)
    plt.show()

def stretch_poi_by_length(values, poi, stretch_by):
    """
    Stretches POI in the values provided to extend the POI region by desired length
    :param values:
    :param poi:
    :param stretch_by:
    :return:
    """
    for poi_name, poi_points in poi.iteritems():
        poi_data = values[poi_points]
        new_len = len(poi_data) + stretch_by
        new_poi = uniform_scaling_to_length(poi_data, new_len)

        # Assume POI regions are contiguous
        before_poi = values[:poi_points[0]]
        after_poi = values[poi_points[-1] + 1:]

        new_data = np.concatenate((before_poi, new_poi, after_poi))
        new_poi = np.asarray(range(poi_points[0], poi_points[0] + new_len), dtype=int)

        return new_data, {poi_name: new_poi}


def read_pois(poi_file, regions, dataset, resolution, ignore_poi_non_overlaps=False, use_strand_information=True):
    print '> Reading points of interest'

    poi = Regions.from_bed(poi_file)
    poi = poi.as_bins_of(regions, resolution=resolution,
                         ignore_non_overlaps=ignore_poi_non_overlaps,
                         account_for_strand_information=use_strand_information)

    dataset.add_points_of_interest(poi, name=os.path.basename(poi_file))
    poi_dataset = dataset.drop_no_pois()

    return poi_dataset

def generate_new_dataset(dataset, desired_size, max_stretch_length=20, reverse_randomly=True, mutation_probability=0.01,
                         mutation_scale=0.2):
    """
    Generates new dataset of desired size by randomly stretching the POI regions of the dataset provided up to
    max_stretch_length

    :param dataset:
    :param desired_size:
    :param max_stretch_length:
    :param reverse_randomly: will reverse some of the regions randomly if set to true
    :param mutation_probability: roughly the percentage of items in the sequence  that will mutate and obtain the value
                                 of some other random index (i.e. 0.01 means that around 1% of items will mutate)
    :param mutation_scale: the scale of mutations to occur. a value of 0.2. means that the value could change from -20% to 20%

    :return:
    """
    points_of_interest = dataset.points_of_interest
    counts = defaultdict(lambda: 0)

    new_dataset = {}
    new_points_of_interest = {}
    for i in range(desired_size):
        item = random.choice(dataset.items)
        data = dataset.ix[item]
        values = data.values
        poi = points_of_interest[item]

        stretch_length = random.randint(0, max_stretch_length)
        new_data, new_poi = stretch_poi_by_length(values, poi, stretch_length)

        new_i = counts[item] + 1
        counts[item] = new_i

        new_ix = '{0}-{1}'.format(item, new_i)

        if reverse_randomly and random.choice([True, False]):
            new_data = reverse_sequence(new_data)
            new_poi[new_poi.keys()[0]] = np.asarray(reverse_sequence((len(new_data) - 1 - new_poi[new_poi.keys()[0]])),
                                                    dtype=int)
            new_ix += 'r'


        if mutation_probability is not None:
            new_data = mutate_sequence(new_data, mutation_probability, mutation_scale)

        new_data = pd.DataFrame(new_data, columns=data.columns)
        new_dataset[new_ix] = new_data
        new_points_of_interest[new_ix] = new_poi


    new_dataset = pd.Panel(new_dataset)
    ad = AlignmentsData(new_dataset)
    ad.points_of_interest = new_points_of_interest

    return ad

def save_dataset(dataset, output_file):
    """
    Will save the dataset to output file that is ready to be passed in DGW as --processed-dataset
    :param dataset:
    :param output_file
    :return:
    """
    f = open(output_file, 'w')
    try:
        pickle.dump(dataset, f)
    except Exception:
        f.close()

def read_regions(regions_filename, resolution):
    regions = Regions.from_bed(regions_filename)
    total_len = len(regions)
    print '> {0} regions of interest read'.format(total_len)

    regions = regions.clip_to_resolution(resolution)

    return regions

if __name__ == '__main__':
    parser = argument_parser()
    args = parser.parse_args()

    print '> Reading regions'
    regions = read_regions(args.regions, args.resolution)
    print '> Reading (full) dataset'
    dataset = read_datasets(args, regions)

    # Convert to log early
    dataset = dataset.to_log_scale()

    print '> Reading POI regions'
    dataset = read_pois(args.points_of_interest, regions, dataset, args.resolution,
                        ignore_poi_non_overlaps=args.ignore_poi_non_overlaps,
                        use_strand_information=args.use_strand_information)

    print '> Picking {0} regions at random'.format(args.number_of_regions)
    random_dataset = get_random_sample(dataset, args.number_of_regions)
    items = random_dataset.items
    print '> Saving regions that were picked at random to {0}'.format(args.output_regions)
    regions.ix[items].to_bed(args.output_regions)

    print '> Generating new dataset'
    new_dataset = generate_new_dataset(random_dataset, args.size,
                                       max_stretch_length=args.max_stretch_length,
                                       reverse_randomly=not args.no_reverse,
                                       mutation_probability=args.mutation_probability,
                                       mutation_scale=args.mutation_scale)

    print '> Saving the newly-generated dataset to {0}'.format(args.output)
    save_dataset(new_dataset, args.output)

    print '> Use --processed-dataset parameter in dgw-worker to perform experments on this dataset'
