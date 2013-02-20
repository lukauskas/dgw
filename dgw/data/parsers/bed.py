import pandas as pd
from dgw.data.containers import Regions

__author__ = 'saulius'
BED_COLUMNS = ['chromosome', 'start', 'end', 'name', 'score',
               'strand', 'thick_start', 'thick_end', 'item_rgb', 'block_count', 'block_sizes', 'block_starts']

def read_bed(bed_file):
    """
    Parses the bed file specified into `pd.DataFrame`
    :param bed_file:
    :return:
    :rtype: `pd.DataFrame`
    """
    regions = pd.read_csv(bed_file, sep="\t", header=None)

    regions.columns = BED_COLUMNS[:len(regions.columns)]
    if 'name' in regions.columns:
        regions = regions.set_index('name')

    return Regions(regions)

def write_bed(regions, file_object):
    import csv

    writer = csv.writer(file_object, delimiter='\t', quotechar='"')
    for ix, data in regions.iterrows():

        row = list(data[['chromosome', 'start', 'end']])
        row.append(ix) # Index == name is the fourth positional argument

        # Score is the fifth
        if 'score' in data:
            row.append(data['score'])

        if 'strand' in data:
            if not 'score' in data:
                row.append(0)  # We do not known the score, so make up a value

            row.append(data['strand'])

        writer.writerow(row)
