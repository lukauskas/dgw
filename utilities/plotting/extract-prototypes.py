import argparse
from dgw.cli import StoreFilenameAction
from dgw.cli.configuration import load_configuration_from_file

__author__ = 'saulius'

def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Use various flags to switch between input types
    parser.add_argument('configuration_file', metavar='dgw_config_file.dgw',
                        help='DGW configuration file outputted by dgw-worker.py',
                        action=StoreFilenameAction)

    parser.add_argument('output_directory', help='Output directory')
    return parser

def main():
    parser = argument_parser()
    args = parser.parse_args()
    configuration = load_configuration_from_file(args.configuration_file)

