__author__ = 'saulius'

import pandas as pd
from dgw.data.containers import Regions

class Genes(Regions):

    def transcription_start_sites(self):
        """
        Returns transcription start site locations for the current set of Genes.
        Pays attention to the strand of the gene and returns the start of the gene for positive strand genes
        and the end of the gene for the negative strand genes.

        :rtype: `pd.Series`
        """
        start_sites = self[self.strand=='+']['start'].append(self[self.strand=='-']['end'])
        start_sites = start_sites.ix[self.index] # reindex the data using original order

        return start_sites

    def regions_around_transcription_start_sites(self, window_width):
        """
        Returns a `Regions` object corresponding to the locations from -window_width to +window_width around the
        transcription start sites for these genes

        :param window_width: the width of the window
        :type window_width: int
        :return: regions around tss
        :rtype: `Regions`
        """

        tss_locations = self.transcription_start_sites()

        # Add window around these
        starts = tss_locations - window_width
        ends   = tss_locations + window_width

        regions_df = pd.DataFrame({'chromosome' : self['chromosome'],
                                'start' : starts,
                                'end' : ends}, index=self.index)

        return Regions(regions_df)


def read_known_genes(known_genes_filename):
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