__author__ = 'saulius'

import pandas as pd
import helpers

def read_known_genes_file(known_genes_filename):

    known_genes = pd.read_csv(known_genes_filename, sep='\t')

    cols = list(known_genes.columns)
    cols = map(lambda x : x.strip('#'), cols)

    known_genes.columns = cols

    return known_genes

def get_tss_peak_df(known_genes, window, resolution=1):
    known_genes = known_genes[['chrom', 'txStart']].drop_duplicates()

    starts = known_genes['txStart'] - window
    ends   = known_genes['txStart'] + window


    peak_df = pd.DataFrame({'chromosome' : known_genes['chrom'],
                            'start' : starts,
                            'end' : ends})
    peak_df = helpers.clip_to_fit_resolution(peak_df, resolution=resolution)

    return peak_df

