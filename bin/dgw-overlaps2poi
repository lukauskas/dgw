#!/usr/bin/env python
import argparse
import sys
from dgw.cli.actions import StoreFilenameAction
from dgw.data.containers import Regions

def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="This application is able to extract regions in poi_filename "
                                                 "that overlap with regions in the input_filename and return them "
                                                 "in DGW-readable format that can then be used to perform DGW matching")

    # Use various flags to switch between input types
    parser.add_argument('input_filename', metavar='main_regions_of_interest.bed',
                        help='Location of the main regions of interest to be parsed', action=StoreFilenameAction)
    parser.add_argument('poi_filename', metavar='points_of_interest.bed',
                        help='Specific points of interest regions that need to be marked as POI in these regions')

    parser.add_argument('-O', '--output', help='Output file', type=argparse.FileType('w'), default=sys.stdout)
    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()

    input_filename = args.input_filename
    poi_filename = args.poi_filename

    input_regions = Regions.from_bed(input_filename)
    poi_regions = Regions.from_bed(poi_filename)

    output_file = args.output
    for ix, region in input_regions.iterrows():
        pois_in_region = poi_regions.contained_within(region)
        if len(pois_in_region) == 0:
            continue
        output_file.write('{0}:{1}\n'.format(ix, pois_in_region.as_printable_list_of_pois()))

    output_file.close()

if __name__ == '__main__':
    main()
