#!/usr/bin/env python
import argparse
import os
import numpy as np
import cPickle as pickle
import dgw
from dgw.cli import Configuration
from dgw.cluster import add_path_data
import dgw.cluster.visualisation


def load_from_pickle(f):
    if isinstance(f, basestring):
        f = open(f, 'rb')
        close = True
    else:
        close = False
    try:
        return pickle.load(f)
    finally:
        if close:
            f.close()

def strict_load(filename):
    try:
        data = load_from_pickle(filename)
        return data
    except Exception, e:
        raise Exception("Failed to read {0!r}, got {1!r}. "
                        "Make sure that you have all files output from DGW in the same directory as the config file"
                         .format(filename, e))


def argument_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('configuration_file', metavar='dgw_config_file.dgw', type=argparse.FileType('r'))
    parser.add_argument('-O', '--output', metavar='output_directory', default="output")

    return parser

def main():
    parser = argument_parser()
    args = parser.parse_args()

    try:
        configuration = Configuration.from_json(args.configuration_file)
    except Exception, e:
        parser.error('Error opening configuration file provided: {0!r}'.format(e))
    finally:
        args.configuration_file.close()

    if not isinstance(configuration, Configuration):
        parser.error('Invalid configuration file provided. Make sure you are specifying the right file')

    configuration.directory = os.path.dirname(args.configuration_file.name)

    if not configuration.linkage_filename:
        parser.error('No linkage filename provided, cannot explore a --blank run')

    if configuration.parsed_regions_filename:
        regions = strict_load(configuration.parsed_regions_filename)
    else:
        regions = None

    dataset = strict_load(configuration.dataset_filename)

    if configuration.linkage_filename:
        try:
            linkage = np.load(configuration.linkage_filename)
        except Exception, e:
            parser.error('Error reading linkage file {0!r}, got {1!r}',format(linkage, e))

        prototypes = strict_load(configuration.prototypes_filename)
        warping_paths = strict_load(configuration.warping_paths_filename)

        print "> Processing data"
        hc = dgw.cluster.analysis.HierarchicalClustering(dataset, regions, linkage_matrix=linkage, prototypes=prototypes,
                                                         dtw_function=configuration.dtw_function,
                                                         prototyping_method=configuration.prototyping_method)
        add_path_data(hc.tree_nodes_list, hc.num_obs, warping_paths)
        configuration_basename = os.path.basename(args.configuration_file.name)
        hcv = dgw.cluster.visualisation.HierarchicalClusteringViewer(hc, output_directory=args.output,
                                                                     configuration_file=configuration_basename)
        print "> Displaying explorer"
        hcv.show()

if __name__ == '__main__':
    main()