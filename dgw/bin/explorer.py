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

    parser.add_argument('--save-only', '-s',
                        help='Do not start the interactive viewer, just save the output.'
                             'Requires output directory and --cut or --n-clusters set',
                        action='store_const', const=True, default=False)
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

    # Borrowed from colorbrewer's Dark2 color palette
    standard_highlight_colours = ["#d95f02", "#e7298a"]
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

            if not poi:
                raise Exception('POI file provided, but no POIs were parsed from {}. Try using dgw-overlaps2poi'.format(poi_file))

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
    have_cut = False
    if args.cut:
        cut_xdata = args.cut
        have_cut = True
    elif args.n_clusters:
        cut_xdata = hc.distance_threshold_for_n_clusters(args.n_clusters)
        have_cut = True

    hcv = dgw.cluster.visualisation.HierarchicalClusteringViewer(hc, output_directory=args.output,
                                                                 configuration_file=configuration_basename,
                                                                 highlight_colours=highlight_colours,
                                                                 cut_xdata=cut_xdata)
    if not args.save_only:
        print "> Displaying explorer"
        hcv.show()
    else:
        if not have_cut:
            raise Exception('Please use specify either the cut distance, or number of clusters when using --save-only')
        output_directory = args.output
        if not output_directory:
            raise Exception('Please specify output directory where the files should be saved')

        if not os.path.isdir(output_directory):
            os.makedirs(output_directory)

        print('> Saving to {}'.format(output_directory))
        hcv.savefig(os.path.join(output_directory, 'clustering.pdf'))

        cluster_previewer = hcv.cluster_previewer()
        cluster_previewer.save_clusters(output_directory)

        print('> Saving summaries')
        cluster_previewer.save_previewer_windows(os.path.join(output_directory, 'summaries'))






if __name__ == '__main__':
    main()