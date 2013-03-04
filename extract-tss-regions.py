#!/usr/bin/env python
import argparse
import sys
from dgw.cli.actions import StoreFilenameAction
from dgw.data.containers import Genes

def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Use various flags to switch between input types
    parser.add_argument('input_filename', help='Location of ENCODE knownGenes file to be parsed', action=StoreFilenameAction)
    parser.add_argument('-w', '--window', help='Window around TSS to consider', type=int, default=2000)
    parser.add_argument('-O', '--output', help='Output file', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-OT', '--tss-output', help='Transcription start sites output', type=argparse.FileType('w'), default=None)

    return parser

def main():
    parser = argument_parser()
    args = parser.parse_args()

    input_filename = args.input_filename

    # Reading genes
    genes = Genes.from_encode_known_genes(input_filename)
    tss_regions = genes.regions_around_transcription_start_sites(args.window)
    tss_regions.to_bed(args.output)

    if args.tss_output:
        tss = genes.regions_around_transcription_start_sites(0)
        tss.data['end'] += 1 # Add one to the end coordinate so regions are of length 1.
        tss.to_bed(args.tss_output)


if __name__ == '__main__':
    main()
