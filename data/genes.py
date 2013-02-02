__author__ = 'saulius'

import pandas as pd

def read_known_genes(known_genes_filename):
    '''
    Reads known_genes file.
    The file can be obtained from the
    ``table browser from ENCODE project <http://encodeproject.org/cgi-bin/hgTables?hgsid=320204017&clade=mammal&org=Human&db=hg19&hgta_group=genes&hgta_track=wgEncodeRegTxn&hgta_table=0&hgta_regionType=genome&position=chr21%3A33%2C031%2C597-33%2C041%2C570&hgta_outputType=wigData&hgta_outFileName=>``

    :param known_genes_filename:
    :return:
    :rtype : pd.DataFrame
    '''

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


    return known_genes

def tss_locations(known_genes):
    '''
    Returns locations of TSS. Either uses gene['start'] or gene['end'] depending on whether the gene is on '+' or '-' strand.

    :param known_genes:
    :return:
    '''

    locations = known_genes[known_genes.strand=='+']['start'].append(known_genes[known_genes.strand=='-']['end'])
    locations = locations.ix[known_genes.index] # reindex the data using original order

    return locations

def get_regions_around_tss(known_genes, window, resolution=1):
    '''
    Returns regions that are around transcription starting sites.

    :param known_genes: parsed known_genes file, see read_known_genes_file
    :param window: window of the number of base pairs to add to region. A window of 2000 would return from -2000 to +2000
                   base pairs around the TSS
    :param resolution: Resolution of data to clip the window around.
    :return:
    '''
    # Get the locations of tss. Note that they coinside with Transcription End Site for peaks on negative strand
    tss_locations = tss_locations(known_genes)

    # Add window around these
    starts = tss_locations - window
    ends   = tss_locations + window

    peak_df = pd.DataFrame({'chromosome' : known_genes['chromosome'],
                            'start' : starts,
                            'end' : ends}, index=known_genes.index)

    peak_df = helpers.clip_to_fit_resolution2(peak_df, resolution=resolution)
    peak_df = peak_df.drop_duplicates()
    return peak_df

def genes_that_may_overlap_region(known_genes, region, overlap_window=2000):
    chr_genes = known_genes[known_genes.chromosome == region['chromosome']]

    transcription_starts = tss_locations(chr_genes)
    overlaps = transcription_starts[transcription_starts.between(region['start'] - overlap_window,
                                                                     region['end'] + overlap_window)]

    overlap_indices = overlaps.index - [region.name] # Remove the region provided from list

    return known_genes.ix[overlap_indices]

def overlapping_gene_counts(known_genes, regions, overlap_window=2000):

    ixs = []
    counts = []

    for ix, region in regions.iterrows():
        ixs.append(ix)
        counts.append(len(genes_that_may_overlap_region(known_genes, region, overlap_window=overlap_window)))

    return pd.Series(counts, index=ixs)

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

    return data