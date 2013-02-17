import pandas as pd
from dgw.data.containers import Regions

__author__ = 'saulius'

def read_bed(bed_file):
    """
    Parses the bed file specified into `pd.DataFrame`
    :param bed_file:
    :return:
    :rtype: `pd.DataFrame`
    """
    regions = pd.read_csv(bed_file, sep="\t", header=None)

    regions.columns = ['chromosome', 'start', 'end', 'name', 'score']
    regions = regions.set_index('name')

    return Regions(regions)