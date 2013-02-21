__author__ = 'saulius'

import pandas as pd
from ..containers import Genes

def read_encode_known_genes(known_genes_filename):
    """
    Reads known_genes file.
    The file can be obtained from the
    ``table browser from ENCODE project <http://encodeproject.org/cgi-bin/hgTables?hgsid=320204017&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=wgEncodeRegTxn&hgta_table=0&hgta_regionType=genome&position=chr21%3A33%2C031%2C597-33%2C041%2C570&hgta_outputType=wigData&hgta_outFileName=>``

    :param known_genes_filename:
    :return:
    :rtype: pd.DataFrame
    """

    known_genes = pd.read_csv(known_genes_filename, sep='\t')

    cols = list(known_genes.columns)
    cols = map(lambda x : x.strip('#'), cols)

    # Rename cols
    new_cols = []
    for c in cols:
        if c == 'chrom':
            c = 'chromosome'
        elif c == 'txStart':
            c = 'start'
        elif c == 'txEnd':
            c = 'end'
        new_cols.append(c)

    known_genes.columns = new_cols
    known_genes = known_genes.set_index('name')

    return Genes(known_genes)

def read_gtf(gtf_location):
    '''
        Reads GTF File
    :param gtf_location:
    :return:
    '''

    def parse_attribute(attribute_str):
        attribute = attribute_str.split(';')
        attribute = map(lambda x: x.strip(), attribute)
        attribute = map(lambda x: x.split(' '), attribute)
        attribute = filter(lambda x : len(x) == 2, attribute)
        attribute = map(lambda x: (x[0], x[1].strip('"')), attribute)

        return dict(attribute)

    data = pd.read_csv(gtf_location, sep='\t', skiprows=5, header=None)
    data.columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

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

    return Genes(data)