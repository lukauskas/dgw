from logging import debug
import pandas as pd
import numpy as np
from dgw.data.parsers.pois import map_to_bins
from dgw.dtw import no_nans_len


class AlignmentsDataIndexer(object):
    """
    A wrapper around `_NDFrameIndexer` that would return `AlignmentsData` objects instead of `pd.Panel` objects
    """
    _ndframe_indexer = None
    _alignments_data = None

    def __init__(self, ndframe_indexer, alignments_data):
        self._ndframe_indexer = ndframe_indexer
        self._alignments_data = alignments_data

    def __getitem__(self, key):
        result = self._ndframe_indexer.__getitem__(key)
        if isinstance(result, pd.Panel):
            data = AlignmentsData(result, self._alignments_data.resolution)
            data.points_of_interest = self._alignments_data.points_of_interest
            return data
        else:
            return result

    def __setitem__(self, key, value):
        self._ndframe_indexer.__setitem__(self, key, value)

class RegionsIndexer(object):
    """
    A wrapper around `_NDFrameIndexer` that would return `Regions` objects instead of `pd.DataFrame` objects
    """
    _ndframe_indexer = None
    regions = None

    def __init__(self, ndframe_indexer, regions):
        self._ndframe_indexer = ndframe_indexer
        self.regions = regions

    def __getitem__(self, key):
        result = self._ndframe_indexer.__getitem__(key)
        if isinstance(result, pd.DataFrame):
            data = self.regions.__class__(result)
            return data
        else:
            return result

    def __setitem__(self, key, value):
        self._ndframe_indexer.__setitem__(self, key, value)

class AlignmentsData(object):

    _data = None
    _poi = None
    _scale = None
    _resolution = None

    def __init__(self, panel, resolution, poi=None, scale='raw'):
        """
        Initialises `AlignmentsData` with a `panel` provided.
        The panel is assumed to have data sets on the minor axis
        See `dgw.data.parsers.read_bam` for how to generate this data.

        :param panel: `pd.Panel` object `AlignmentsData` will be initialised on, or `pd.DataFrame` that will be converted
                     to Panel
        :param resolution: resolution of data
        :param poi: points of interest
        :param scale: the scale of data

        :return:
        """
        if isinstance(panel, pd.DataFrame):
            # Create a panel from the DataFrame by giving it a generic name and making sure it is on the minor axis
            self._data = pd.Panel(set(['Dataset 1'])).transpose(1, 2, 0)
        elif isinstance(panel, pd.Panel):
            self._data = panel
        else:
            raise Exception('Invalid type of data provided for AlignmentsData: {0!r}, expected pd.Panel or pd.Dataframe'
                            .format(type(panel)))

        self._scale = scale

        self.points_of_interest = poi
        self._resolution = resolution


    def reset_poi(self):
        self._poi = {}

    @property
    def data(self):
        return self._data

    @property
    def points_of_interest(self):
        return self._poi

    @property
    def resolution(self):
        return self._resolution

    @points_of_interest.setter
    def points_of_interest(self, value):
        if value is None:
            self._poi = {}
        else:
            self._poi = dict([(ix, value) for ix, value in value.iteritems() if ix in self.items])


    def add_points_of_interest(self, binned_points_of_interest_regions, name):
        for ix, value in binned_points_of_interest_regions.iteritems():
            if ix not in self.items:
                continue

            try:
                self._poi[ix][name] = value
            except KeyError:
                self._poi[ix] = {name: value}


    def drop_no_pois(self):
        common_index = self.items & self.points_of_interest

        return self.ix[common_index]

    #-- Functions that simulate pd.Panel behaviour -------------------------------------------------------------------

    def mean(self, axis='items', skipna=True):
        # Override axis parameter in the pd.Panel mean function
        return self._data.mean(axis=axis, skipna=skipna)

    @property
    def values(self):
        return self.data.values


    @property
    def number_of_datasets(self):
        return len(self.dataset_axis)

    @property
    def number_of_items(self):
        return len(self.items)

    @property
    def number_of_columns(self):
        return len(self.major_axis)

    @property
    def lengths(self):
        a = []
        for _, row in self.data.iteritems():
            a.append(no_nans_len(row.values))

        return pd.Series(a, index=self.items)


    @property
    def dataset_axis(self):
        return self.data.minor_axis

    def dataset_xs(self, *args, **kwargs):
        return self.data.minor_xs(*args, **kwargs)

    @property
    def items(self):
        return self.data.items

    @property
    def ix(self):
        return AlignmentsDataIndexer(self.data.ix, self)

    def __getitem__(self, item):
        return self.data.__getitem__(item)

    def head(self, n=5):
        return self.ix[:n]

    @property
    def major_axis(self):
        return self.data.major_axis


    #-- Additional transformations not visible in default pd.Panel  ----------------------------------------
    def __len__(self):
        return self.number_of_items

    def to_log_scale(self):
        if self._scale == 'log':
            return self

        new_data = (self.data + 2).apply(np.log)  # Adding +2 so we have no zeros in log output
        ad = AlignmentsData(new_data, self.resolution, scale='log')
        ad.points_of_interest = self.points_of_interest
        return ad

    def normalise_bin_heights(self):
        data = {}
        for ix, data_row in self.data.iteritems():
            data[ix] = data_row / data_row.max()

        data = pd.Panel(data)
        return self.__class__(data, self.resolution, poi=self.points_of_interest)

    def plot_heatmap(self, *args, **kwargs):
        """
        Plots heatmap of the data stored in the panel.

        :param args: args to be passed into `visualisation.heatmap.plot`
        :param kwargs: kwargs to be passed into `visualisation.heatmap.plot`
        :return:
        """
        import visualisation.heatmap
        return visualisation.heatmap.plot(self, *args, **kwargs)


    def __repr__(self):
        return '{0} containing\n{1!r}'.format(self.__class__, self.data)

    def __array__(self, *args, **kwargs):
        return self.data.__array__(*args, **kwargs)

class Regions(object):
    REQUIRED_COLUMNS = frozenset(['chromosome', 'start', 'end'])

    _data = None
    def __init__(self, data, *args, **kwargs):

        if isinstance(data, Regions):
            self._data = data.data
        else:
            data = pd.DataFrame(data, *args, **kwargs)
            # Verify that all required columns are in the DF
            for column in self.REQUIRED_COLUMNS:
                if column not in data.columns:
                    raise ValueError('No such column {0!r} in provided DataFrame'.format(column))

            self._data = data

        self._data = self._data.drop_duplicates() # TODO: somehow join the indices of the dropped data

    @property
    def data(self):
        return self._data

    def __repr__(self):
        return '{0} containing \n{1!r}'.format(self.__class__, self.data)


    # --- Initialisation ----------------------------------------------------------------------------------------------
    @classmethod
    def from_bed(cls, bed_file):
        from dgw.data.parsers import read_bed
        return cls(read_bed(bed_file))

    def to_bed(self, bed_file, **track_kwargs):
        from dgw.data.parsers import write_bed
        return write_bed(self, bed_file, **track_kwargs)

    # --- Functions that provide direct access to the DataFrame behind all this ----------------------------------------
    def __getitem__(self, item):
        result = self.data[item]
        if isinstance(result, pd.DataFrame):
            try:
                # Try returining it as the same class
                return self.__class__(result)
            except ValueError:
                # If not valid, return as DataFrame
                return result
        return result

    def has_strand_data(self):
        return 'strand' in self.columns

    def infer_strand_from_whether_the_region_was_reversed_or_not(self, reversion_status_dictionary):
        """
        Infers the strand of a region by checking whether the region was reversed or not.

        :param reversion_status_dictionary: Dictionary of boolean true/false values, true indicated the region
                                            was reversed by DTW
        :return:
        """
        data = self.data.copy()

        strand_dict = {}
        for key, reversed in reversion_status_dictionary.iteritems():
            if key not in data.index:
                continue
            if reversed is None:
                strand_dict[key] = None
            elif reversed:
                strand_dict[key] = '-'
            else:
                strand_dict[key] = '+'

        strand_series = pd.Series(strand_dict)
        data['strand'] = strand_series
        return self.__class__(data)

    def iterrows(self):
        return self.data.iterrows()

    def head(self, *args, **kwargs):
        return self.__class__(self.data.head(*args, **kwargs))

    @property
    def index(self):
        return self.data.index

    def join(self, *args, **kwargs):
        new_data = self.data.join(*args, **kwargs)
        return self.__class__(new_data)

    def append(self, *args, **kwargs):
        new_data = self.data.append(*args, **kwargs)
        return self.__class__(new_data)

    @property
    def ix(self):
        return RegionsIndexer(self.data.ix, self)

    @property
    def columns(self):
        return self.data.columns

    def __getattr__(self, name):
        """
        Emulate the behaviour in `pd.DataFrame` to return one of the columns as an attribute.
        :param name:
        :return:
        """
        # Partially stolen from pandas implementation.
        if name in self.data.columns:
            return self.data[name]

        raise AttributeError("{0!r} has no attribute {1!r}".format(type(self).__name__, name))

    def __len__(self):
        return self.data.__len__()

    # --- Functions special to Regions ---------------------------------------------------------------------------------
    @property
    def lengths(self):
        """
        Returns a `pd.Series` containing the lengths of all the regions contained
        :rtype: `pd.Series`
        """
        series = self['end'] - self['start']
        series.name = 'length'
        return series

    def contained_within(self, other_region):
        """
        Returns all regions that are contained within other region.
        Returns only those regions that are fully contained within the query region.

        :param other_region: The region that regions will be checked to be inside
        :type other_region: `pd.Series`
        :return:
        """
        return self[(self.chromosome == other_region['chromosome'])
                    & (self.start >= other_region['start'])
                    & (self.end <= other_region['end'])]

    def as_printable_list_of_pois(self):

        printable_list = ""
        for ix, row in self.iterrows():
            if printable_list:
                printable_list += ','
            printable_list += ','.join(map(str, range(row['start'], row['end'])))

        return printable_list

    def regions_not_in_dataset(self, dataset):
        """
        Returns all regions that are not in the dataset provided.

        :param dataset:
        :type dataset: AlignmentsData
        :return:
        """
        missing_indices = self.index[~self.index.isin(dataset.items)]
        return Regions(self.data.ix[missing_indices])

    def as_bins_of(self, other_regions, resolution=1, ignore_non_overlaps=False, account_for_strand_information=False):
        """
        Returns the current regions as bins of some other set of regions provided
        :param other_regions: regions to match self to
        :param resolution: resolution at which to do so
        :param ignore_non_overlaps: if set to false, the parser will raise a ValueError if self does not overlap with other regions
        :param account_for_strand_information: if set to true, the parser will account for the antisense regions and
            return bins from the end, rather than front
        :return:
        """
        other_regions = other_regions.clip_to_resolution(resolution)

        bins = {}
        for ix, data in self.iterrows():
            current_start = data['start']
            current_end = data['end']
            current_chromosome = data['chromosome']

            try:
                other = other_regions.ix[ix]
            except KeyError:
                continue

            other_chromosome = other['chromosome']

            if current_chromosome != other_chromosome:
                if ignore_non_overlaps:
                    continue
                else:
                    raise ValueError('Points of interest do not overlap with regions of interest. Failed ix:{0!r}'.format(ix))

            bins[ix] = map_to_bins(range(current_start, current_end), other, resolution=resolution,
                                   ignore_non_overlaps=ignore_non_overlaps,
                                   account_for_strand_information=account_for_strand_information)

        return bins



    def clip_to_resolution(self, resolution):
        """
        Adjusts the region boundaries to fit the resolution by extending the regions boundaries as to length of the
        regions is divisible from resolution.

        Returns a new `Regions` object rather than clipping this one in place

        Please note that the new `Regions` object this function returns will not include any other fields but
        the ones listed in `Regions.REQUIRED_COLUMNS`.

        :param resolution: the resolution to clip
        :type resolution: int
        :return: a new `Regions` object
        :rtype: Regions
        """

        resolution = int(resolution)

        if resolution == 1:
            return self
        elif resolution <= 0:
            raise ValueError('Resolution should be > 0 ({0!r} given)'.format(resolution))

        lengths = self.lengths
        new_regions_data = []

        for ix, row in self.iterrows():
            length = lengths.ix[ix]
            remainder = length % resolution

            if remainder == 0:
                offset_needed = 0
            else:
                offset_needed = resolution - remainder

            add_left = offset_needed / 2
            add_right = offset_needed / 2 + offset_needed % 2

            row['start'] -= add_left

            if row['start'] < 0:
                # Check if we accidentally went sub zero
                add_right += -row['start']
                row['start'] = 0

            row['end'] += add_right

            new_regions_data.append(row)

        df = pd.DataFrame(new_regions_data, index=self.index)
        return Regions(df)

    def bins_to_intervals(self, resolution):
        """
        Returns genomic location intervals in the form of (start, end) for each bin of the region at resolution
        provided, for each region
        :param resolution:
        :return:
        """
        regions = self.clip_to_resolution(resolution)

        intervals = {}
        for ix, row in regions.iterrows():
            start = row['start']
            end = row['end']

            intervals[ix] = [(i, i+resolution) for i in np.arange(start, end, resolution, dtype=int)]

        return intervals



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
        starts[starts < 0] = 0
        ends = tss_locations + window_width

        regions_df = pd.DataFrame({'chromosome' : self.chromosome,
                                   'start' : starts,
                                   'end' : ends,
                                   'strand' : self.strand}, index=self.index)

        return Regions(regions_df)

    def regions_around_splicing_site(self, splicing_site, window_width):
        """
        Return a `Regions` object with windows around first splicing site
        :param splicing_site: number of the splicing site to get (0 based - "0" for first splicing site)
        :param window_width:
        :return:
        """
        exon = self.get_exon_regions(splicing_site)

        exon_genes = self.ix[exon.index]

        ss = exon[exon_genes.strand == '+']['end'].append(exon[exon_genes.strand == '-']['start'])
        ss = ss.ix[self.index]

        starts = ss - window_width
        ends = ss + window_width

        regions_df = pd.DataFrame({'chromosome': exon_genes['chromosome'],
                                  'start': starts, 'end': ends, 'strand' : exon_genes.strand},
                                  index=exon_genes.index)

        return Regions(regions_df)

    def get_exon_regions(self, exon_number, account_for_antisense=True):
        """
        Returns exon regions for particular exon_number provided (0-based)
        :param exon_number: exon number - 0 based. I.e. use 0 to get regions of first exon
        :param account_for_antisense: if set to true, it will automatically account for antisense genes and return correct exon
        :return:
        """
        if 'exonStarts' not in self.columns or 'exonEnds' not in self.columns:
            raise Exception('No exon data available, sorry')

        sub_df = self[['chromosome', 'exonStarts', 'exonEnds', 'strand']]
        new_df = []
        for ix, row in sub_df.iterrows():
            exon_starts = row['exonStarts'].strip(',').split(',')
            exon_ends = row['exonEnds'].strip(',').split(',')
            chromosome = row['chromosome']
            strand = row['strand']

            if account_for_antisense and strand == '-':
                exon_i = len(exon_starts) - 1 - exon_number
            else:
                exon_i = exon_number

            try:
                start = int(exon_starts[exon_i])
                end = int(exon_ends[exon_i])
            except IndexError:
                start = np.nan
                end = np.nan

            new_df.append({'chromosome': chromosome, 'start': start, 'end': end, 'strand': strand})
        new_df = pd.DataFrame(new_df, index=self.index).dropna()
        return Regions(new_df)

    @classmethod
    def from_encode_known_genes(cls, encode_known_genes_filename):
        from dgw.data.parsers import read_encode_known_genes
        return read_encode_known_genes(encode_known_genes_filename)

    @classmethod
    def from_gtf(cls, gtf_filename):
        from dgw.data.parsers import read_gtf
        return read_gtf(gtf_filename)
