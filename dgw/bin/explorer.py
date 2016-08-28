#!/usr/bin/env python
import argparse
import os
import dgw
from dgw.cli import StoreUniqueFilenameAction
from dgw.cli.configuration import load_configuration_from_file
import dgw.cluster.visualisation
from dgw.data.containers import Regions
from dgw.data.parsers.pois import from_simple


def argument_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('configuration_file', metavar='dgw_config_file.dgw', type=argparse.FileType('r'))
    parser.add_argument('-O', '--output', metavar='output_directory', default="output")
    poi_group = parser.add_mutually_exclusive_group()
    poi_group.add_argument('-poi', '--points-of-interest', metavar='poi.bed', nargs='+', action=StoreUniqueFilenameAction,
                        help='Points of interest. Note if the dataset has points of interest associated with it, '
                             'they will be overwritten')

    poi_group.add_argument('--no-poi', action='store_const', const=True, default=False,
                        help='Do not display any POI anotations that are bundled with DGW-worker result')
    parser.add_argument('--ignore-poi-non-overlaps', default=False, action='store_const', const=True,
                        help='If set to true, DGW will silently ignore -pois that do not overlap with the regions')

    cut_group = parser.add_mutually_exclusive_group()
    cut_group.add_argument('--cut', '-c', type=float, default=None, help='Cut threshold to initialise DGW to.')
    cut_group.add_argument('--n-clusters', '-nc', type=int,
                           default=None,
                           help='Number of clusters to initialise the dendrogram cut to.')


    parser.add_argument('-v', '--verbose', action='store_const', const=True, default=False)

    return parser

def main():
    parser = argument_parser()
    args = parser.parse_args()

    if args.verbose:
        import logging
        logging.root.setLevel(logging.DEBUG)

    configuration = load_configuration_from_file(args.configuration_file)
    args.configuration_file.close()

    if configuration.blank:
        parser.error('Cannot explore a --blank run of DGW')

    regions = configuration.load_regions()
    dataset = configuration.load_dataset()

    standard_highlight_colours = ["w", "k"]
    highlight_colours = {}
    if args.points_of_interest:
        dataset.reset_poi()
        for i, poi_file in enumerate(args.points_of_interest):
            print '> Reading points of interest from {0!r}'.format(poi_file)

            try:
                poi = from_simple(poi_file, regions, resolution=configuration.resolution,
                                  account_for_strand_information=configuration.use_strand_information)
            except ValueError:

                poi = Regions.from_bed(poi_file)
                poi = poi.as_bins_of(regions, resolution=configuration.resolution,
                                     ignore_non_overlaps=args.ignore_poi_non_overlaps,
                                     account_for_strand_information=configuration.use_strand_information)

            poi_filename = os.path.basename(poi_file)
            dataset.add_points_of_interest(poi, name=poi_filename)
            try:
                highlight_colours[poi_filename] = standard_highlight_colours.pop()
            except IndexError:
                raise Exception("Sorry, only up to {0} POI regions are supported".format(len(standard_highlight_colours)))
    else:
        if args.no_poi:
            dataset.reset_poi()
            
        if dataset.points_of_interest:
            highlight_colours[dataset.points_of_interest.values()[0].keys()[0]] = standard_highlight_colours.pop()


    hc = configuration.create_hierarchical_clustering_object(regions=regions, dataset=dataset)
    configuration_basename = os.path.basename(args.configuration_file.name)

    cut_xdata = 0
    if args.cut:
        cut_xdata = args.cut
    elif args.n_clusters:
        cut_xdata = hc.distance_threshold_for_n_clusters(args.n_clusters)

    hcv = dgw.cluster.visualisation.HierarchicalClusteringViewer(hc, output_directory=args.output,
                                                                 configuration_file=configuration_basename,
                                                                 highlight_colours=highlight_colours,
                                                                 cut_xdata=cut_xdata)

    print "> Displaying explorer"
    hcv.show()

if __name__ == '__main__':
    main()