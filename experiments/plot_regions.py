import argparse
import os
from dgw import read_bam, Regions
from dgw.cli import StoreFilenameAction, StoreUniqueFilenameAction
import matplotlib.pyplot as plt

__author__ = 'saulius'

def read_datasets(args, regions):
    data_filters = []

    return read_bam(args.datasets, regions,
                    resolution=args.resolution,
                    extend_to=args.extend_to,
                    data_filters=data_filters,
                    output_removed_indices=False,
                    reverse_negative_strand_regions=args.use_strand_information)

def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--regions', metavar='regions_of_interest.bed', action=StoreFilenameAction,
                             help='A BED file listing genome regions that should be plotted', required=False)

    parser.add_argument('-d', '--datasets', metavar='dataset.bam', nargs='+', action=StoreUniqueFilenameAction,
                             help='One or more datasets to be plotted. Must be BAM files.')

    parser.add_argument('-res', '--resolution', help='Read resolution', type=int, default=50)
    parser.add_argument('-ext', '--extend_to', help='Extend reads to specified length', type=int, default=200)

    parser.add_argument('-O', '--output-directory', help='Output directory to save regions images to', required=True)

    parser.add_argument('--use-strand-information', const=True, action="store_const", default=False,
                                      help='Will use the strand information provided in the BED file, if set to true.')
    return parser

def read_regions(regions_filename, resolution):
    regions = Regions.from_bed(regions_filename)
    total_len = len(regions)
    print '> {0} regions of interest read'.format(total_len)

    regions = regions.clip_to_resolution(resolution)

    return regions

def main():
    parser = argument_parser()
    args = parser.parse_args()

    print '> Reading regions'
    regions = read_regions(args.regions, args.resolution)
    print '> Loading dataset'
    dataset = read_datasets(args, regions)
    dataset = dataset.to_log_scale()

    output_dir = args.output_directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print '> Saving items'
    for ix, item in dataset.data.iteritems():
        item.plot(legend=False)
        filename = os.path.join(output_dir, '{0}.pdf'.format(ix.replace('.', '_')))
        print 'Saving region {0!r} to filename {1!r}'.format(ix, filename)
        plt.title(ix)
        plt.savefig(filename)
        plt.close('all')

    print '> Done'

if __name__ == '__main__':
   main()