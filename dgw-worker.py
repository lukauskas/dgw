"""
A worker module for Dynamic Genome Warping.

Processes the data into an intermediate representation that could then be analysed.
Designed to be run on a multi-core machine with a lot of memory, e.g. a supercomputer.

"""
import argparse
import os

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


def _argument_parser():

    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--regions', metavar='regions_of_interest.bed', action=StoreFilenameAction,
                        help='A BED file listing genome regions that will be processed', required=True)
    parser.add_argument('-d', '--datasets', metavar='dataset.bam', nargs='+', action=StoreUniqueFilenameAction,
                        help='One or more datasets to be analysed using DGW. Must be BAM files.', required=True)

    parser.add_argument('-p', '--prefix', help='Prefix of the output files generated '
                                               '(defaults to \'dgw\')')

    parser.add_argument('--truncate_regions', metavar='X', type=int, help='Only use first X rows of regions '
                                                                           'rather than the full dataset')

    return parser


def main():
    parser = _argument_parser()
    args = parser.parse_args()
    if args.prefix is None:
        args.prefix = 'dgw'


    print 'Regions: {0!r}'.format(args.regions)
    print 'Datasets: {0!r}'.format(args.datasets)
    print 'Prefix: {0!r}'.format(args.prefix)
if __name__ == '__main__':
    main()