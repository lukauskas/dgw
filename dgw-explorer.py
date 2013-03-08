#!/usr/bin/env python
import argparse
from collections import defaultdict
import os
import numpy as np
import cPickle as pickle
import dgw
from dgw.cli import Configuration, StoreFilenameAction, StoreUniqueFilenameAction
from dgw.cli.configuration import load_configuration_from_file
from dgw.cluster import add_path_data
import dgw.cluster.visualisation
from dgw.data.containers import Regions
import logging




def argument_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('configuration_file', metavar='dgw_config_file.dgw', type=argparse.FileType('r'))
    parser.add_argument('-O', '--output', metavar='output_directory', default="output")
    parser.add_argument('-poi', '--points-of-interest', metavar='poi.bed', nargs='+', action=StoreUniqueFilenameAction,
                        help='Points of interest. Note if the dataset has points of interest associated with it, '
                             'they will be overwritten')

    parser.add_argument('--ignore-poi-non-overlaps', default=False, action='store_const', const=True,
                        help='If set to true, DGW will silently ignore -pois that do not overlap with the regions')

    parser.add_argument('-v', '--verbose', action='store_const', const=True, default=False)

    return parser

def main():
    parser = argument_parser()
    args = parser.parse_args()

    if args.verbose:
        import logging
        logging.root.setLevel(logging.DEBUG)

    configuration = load_configuration_from_file(args.configuration_file)
    configuration.directory = os.path.dirname(args.configuration_file.name)

    if configuration.blank:
        parser.error('Cannot explore a --blank run of DGW')

    regions = configuration.load_regions()
    dataset = configuration.load_dataset()

    if args.points_of_interest:
        combined_poi = defaultdict(lambda: {})

        for i, poi_file in enumerate(args.points_of_interest):
            print '> Reading points of interest from {0!r}'.format(poi_file)
            poi = Regions.from_bed(poi_file)
            poi = poi.as_bins_of(regions, resolution=configuration.resolution,
                                 ignore_non_overlaps=args.ignore_poi_non_overlaps)

            for key, value in poi.iteritems():
                combined_poi[key][i] = value

        dataset.points_of_interest = combined_poi

    hc = configuration.create_hierarchical_clustering_object(regions=regions, dataset=dataset)
    configuration_basename = os.path.basename(args.configuration_file.name)
    hcv = dgw.cluster.visualisation.HierarchicalClusteringViewer(hc, output_directory=args.output,
                                                                 configuration_file=configuration_basename)
    print "> Displaying explorer"
    hcv.show()

if __name__ == '__main__':
    main()