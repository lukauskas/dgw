#!/usr/bin/env python
import argparse
import os
from dgw.cli import StoreFilenameAction
from dgw.cli.configuration import load_configuration_from_file

import matplotlib.pyplot as plt

def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Use various flags to switch between input types
    parser.add_argument('configuration_file', metavar='dgw_config_file.dgw',
                        help='DGW configuration file outputted by dgw-worker.py',
                        action=StoreFilenameAction)

    parser.add_argument('output_directory', help='Output directory')
    return parser

def graph_from_node_list_to_file(node_list, file_object):
    file_object.write('graph dgw_clustering {\n')
    file_object.write('node [style="filled"];')

    for node in node_list:
        file_object.write('"{0}" [label="{0} - {1}", image="{0}.png"];\n'.format(node.id, node.dist))
        if not node.is_leaf():
            file_object.write('"{0}" -- "{1}";\n'.format(node.id, node.get_left().id))
            file_object.write('"{0}" -- "{1}";\n'.format(node.id, node.get_right().id))
    file_object.write('}\n')

def prototype_images_from_node_list(node_list, output_directory):
    print '> Saving images to {0!r}'.format(output_directory)
    ndims = node_list[0].prototype.shape[1]
    plt.figure(figsize=(3 * ndims, 2))
    for node in node_list:
        filename = '{0}.png'.format(node.id)
        full_filename = os.path.join(output_directory, filename)
        print '> Saving {0}'.format(full_filename)
        prototype_T = node.prototype.values.T
        for i in xrange(ndims):
            plt.subplot(1, ndims, i+1)
            plt.plot(prototype_T[i])

        plt.savefig(full_filename)
        plt.clf()

def main():
    parser = argument_parser()
    args = parser.parse_args()
    configuration = load_configuration_from_file(args.configuration_file)

    hc = configuration.create_hierarchical_clustering_object()
    tree_nodes = hc.tree_nodes_list

    output_directory = args.output_directory
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    prototype_images_from_node_list(tree_nodes, output_directory)

    print '> Creating a dot graph'

    f = open(os.path.join(output_directory, 'graph.dot'), 'w')
    graph_from_node_list_to_file(tree_nodes, f)
    f.close()

    print '> Done'

if __name__ == '__main__':
    main()

