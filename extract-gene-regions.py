#!/usr/bin/env python
import argparse
import sys
from dgw.cli.actions import StoreFilenameAction
from dgw.data.containers import Genes

def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Use various flags to switch between input types
    parser.add_argument('input_filename', help='Location of ENCODE knownGenes file to be parsed', action=StoreFilenameAction)
    parser.add_argument('-O', '--output', help='Output file', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-e', '--exon', metavar='X', help='Exon number to output (0-based)', default=0, type=int)
    parser.add_argument('-OE', '--output-exon', metavar='exon_output_filename.bed',
                        help='Filename where exon regions will be saved', type=argparse.FileType('w'))
    return parser

def main():
    parser = argument_parser()
    args = parser.parse_args()

    input_filename = args.input_filename

    # Reading genes
    genes = Genes.from_encode_known_genes(input_filename)
    genes.to_bed(args.output)

    if args.output_exon:
        exons = genes.get_exon_regions(args.exon)
        exons.to_bed(args.output_exon)
        args.output_exon.close()

if __name__ == '__main__':
    main()
