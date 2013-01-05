__author__ = 'saulius'

import pandas as pd
import helpers

def read_known_genes_file(known_genes_filename):

    known_genes = pd.read_csv(known_genes_filename, sep='\t')

    cols = list(known_genes.columns)
    cols = map(lambda x : x.strip('#'), cols)

    known_genes.columns = cols
    known_genes = known_genes.set_index('name')

    return known_genes

def get_tss_peak_df(known_genes, window, resolution=1):
    known_genes = known_genes[['chrom', 'txStart']].drop_duplicates()

    starts = known_genes['txStart'] - window
    ends   = known_genes['txStart'] + window


    peak_df = pd.DataFrame({'chromosome' : known_genes['chrom'],
                            'start' : starts,
                            'end' : ends}, index=known_genes.index)
    peak_df = helpers.clip_to_fit_resolution2(peak_df, resolution=resolution)
    peak_df = peak_df.drop_duplicates()
    return peak_df




def read_gtf(gtf_location):

    def parse_attribute(attribute_str):
        attribute = attribute_str.split(';')
        attribute = map(lambda x: x.strip(), attribute)
        attribute = map(lambda x: x.split(' '), attribute)
        attribute = filter(lambda x : len(x) == 2, attribute)
        attribute = map(lambda x: (x[0], x[1].strip('"')), attribute)

        return dict(attribute)

    data = pd.read_csv(gtf_location, sep='\t', skiprows=5, header=None)
    data.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    # Leave only data about genes
    data = data[data.feature == 'gene']

    # Parse attributes
    parsed_attributes = map(parse_attribute, data['attribute'].values)
    parsed_attributes = pd.DataFrame(parsed_attributes, index=data.index)

    # Join everything into one dataframe
    data = data.join(parsed_attributes)

    # Drop the "attribute" and "feature" columns
    data = data[data.columns - ['attribute', 'feature']]
    data = data.set_index('gene_id')

    # Reorder columns so first_columns are in the beginning of DF
    first_columns = ['seqname', 'start', 'end', 'strand', 'gene_name', 'gene_status', 'gene_type', 'score', 'level']
    data = pd.DataFrame(data, columns=first_columns + list(data.columns - first_columns))

    return data