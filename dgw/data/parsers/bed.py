from StringIO import StringIO
from logging import debug
import pandas as pd
from ..containers import Regions

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
    f = open(bed_file, 'r')
    try:
        s = StringIO()
        # Filter out all lines that do not start with "chr" as BED files are allowed to contain some junk
        for line in f:
            if line.startswith('chr'):
                s.write(line)
        s.seek(0)
        regions = pd.read_csv(s, sep="\t", header=None)
    finally:
        f.close()
        s.close()
    regions.columns = BED_COLUMNS[:len(regions.columns)]

    if len(regions.name) != len(regions.name.drop_duplicates()):
        raise Exception('Input BED file {0!r} contains duplicate values in name column. '
                        'Please ensure the names of the regions are unique'.format(bed_file))

    if 'name' in regions.columns:
        regions = regions.set_index('name')


    return Regions(regions)

def write_bed(regions, writable_file, **track_kwargs):
    import csv

    need_closing = False
    if not isinstance(writable_file, file):
        writable_file = open(writable_file, 'w')
        need_closing = True

    if track_kwargs:
        kwarg_strings = set(["{0}=\"{1}\"".format(key, value) for key,value in track_kwargs.iteritems()])
        track_header = "track {0}\n".format(' '.join(kwarg_strings))
        writable_file.write(track_header)

    writer = csv.writer(writable_file, delimiter='\t', quotechar='"')
    for ix, data in regions.iterrows():

        row = [data['chromosome'], int(data['start']), int(data['end']), ix]

        # Score is the fifth
        if 'score' in data:
            row.append(data['score'])

        if 'strand' in data:
            if not 'score' in data:
                row.append(0)  # We do not known the score, so make up a value

            row.append(data['strand'])

        writer.writerow(row)

    if need_closing:
        writable_file.close()
