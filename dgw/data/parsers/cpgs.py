from ..containers import Regions

__author__ = 'saulius'
import pandas as pd

def read_cpgs(cpgs_filename):
    cpgs = pd.read_csv(cpgs_filename, sep='\t')

    cols = list(cpgs.columns)
    cols = map(lambda x : x.strip('#'), cols)

    # Rename cols
    new_cols = []
    for c in cols:
        if c == 'chrom':
            c = 'chromosome'
        elif c == 'chromStart':
            c = 'start'
        elif c == 'chromEnd':
            c = 'end'
        new_cols.append(c)

    cpgs.columns = new_cols

    return Regions(cpgs)