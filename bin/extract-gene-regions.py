#!/usr/bin/env python
import argparse
import sys
from dgw.cli.actions import StoreFilenameAction
from dgw.data.containers import Genes

def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Use various flags to switch between input types
    parser.add_argument('input_filename', help='Location of ENCODE knownGenes file to be parsed', action=StoreFilenameAction)
    parser.add_argument('output_filename', help='Output file', type=argparse.FileType('w'), default=sys.stdout)

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('-g', '--gene', help='Output whole gene', action='store_const', const=True, default=False)
    group.add_argument('-e', '--exon', metavar='EXON_NUMBER',
                       help='Output the exon EXON_NUMBER (0-based - use "0" for first exon) of each gene',
                       default=None, type=int)
    group.add_argument('-s', '--splicing-site', metavar='EXON_NUMBER',
                       help='Output the splicing site of the exon EXON_NUMBER (0-based again)',
                       default=None, type=int)

    return parser

def main():
    parser = argument_parser()
    args = parser.parse_args()

    input_filename = args.input_filename

    # Reading genes
    genes = Genes.from_encode_known_genes(input_filename)
    output_filename = args.output_filename.name
    if args.gene:
        print '> Saving gene regions to {0}'.format(output_filename)
        genes.to_bed(args.output_filename)
    elif args.exon >= 0:
        exon_regions = genes.get_exon_regions(args.exon)
        print '> Saving all exons {0} to {1}'.format(args.exon, output_filename)
        exon_regions.to_bed(args.output_filename)
    elif args.splicing_site >= 0:
        # We are using hacky way of getting these splicing site locations
        splicing_sites = genes.regions_around_splicing_site(args.splicing_site, 0)
        splicing_sites.data.end += 1  # Window of 0 makes lengths of all regions to be equal to 0, this will make them equal to 1
        print '> Saving all the splicing sites to {0}'.format(output_filename)
        splicing_sites.to_bed(args.output_filename)
    else:
        raise Exception('Invalid arguments provided')

if __name__ == '__main__':
    main()
