import pandas as pd
import matplotlib.pyplot as plt

import dgw.data.visualisation.heatmap as heatmap

class AlignmentsData(object):

    _data = None
    def __init__(self, panel):
        """
        Initialises `AlignmentsData` with a `panel` provided.
        The panel is assumed to have data sets on the minor axis
        See `dgw.data.parsers.read_bam` for how to generate this data.

        :param panel: `pd.Panel` object `AlignmentsData` will be initialised on, or `pd.DataFrame` that will be converted
                     to Panel
        :return:
        """
        if isinstance(panel, pd.DataFrame):
            # Create a panel from the DataFrame by giving it a generic name and making sure it is on the minor axis
            self._data = pd.Panel({'Dataset 1'}).transpose(1, 2, 0)
        elif isinstance(panel, pd.Panel):
            self._data = panel
        else:
            raise Exception('Invalid type of data provided for AlignmentsData: {0!r}, expected pd.Panel or pd.Dataframe'
                            .format(type(panel)))

    @classmethod
    def from_bam(cls, bam_files):
        from dgw.data.parsers import read_bam
        return cls(read_bam(bam_files))

    @property
    def data(self):
        return self._data

    def mean(self, axis='items', skipna=True):
        # Override axis parameter in the pd.Panel mean function
        return self._data.mean(axis=axis, skipna=skipna)

    @property
    def number_of_datasets(self):
        return len(self.dataset_axis)

    @property
    def dataset_axis(self):
        return self.data.minor_axis

    def plot_heatmap(self, *args, **kwargs):
        """
        Plots heatmap of the data stored in the panel.

        :param args: args to be passed into `data.visualisation.heatmap.plot`
        :param kwargs: kwargs to be passed into `data.visualisation.heatmap.plot`
        :return:
        """
        number_of_datasets = self.number_of_datasets

        for i, title in enumerate(self.dataset_axis):
            if number_of_datasets > 1:
                plt.subplot(1, number_of_datasets, i+1) # TODO: consider doing sublot with multiple lines

            data_to_plot = self.dataset_axis(title, copy=False).T
            heatmap.plot(data_to_plot, *args, **kwargs)
            plt.title(title)

    def __repr__(self):
        return '{0} containing\n{1!r}'.format(self.__class__, self.data)

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

    # --- Functions that provide direct access to the DataFrame behind all this ----------------------------------------
    def __getitem__(self, item):
        result = self.data[item]
        if isinstance(result, pd.DataFrame):
            return self.__class__(result)
        return result

    def iterrows(self):
        return self.data.iterrows()

    def head(self, *args, **kwargs):
        return self.__class__(self.data.head(*args, **kwargs))

    def index(self, *args, **kwargs):
        return self.data.index(*args, **kwargs)

    def join(self, *args, **kwargs):
        new_data = self.data.join(*args, **kwargs)
        return self.__class__(new_data)

    def append(self, *args, **kwargs):
        new_data = self.data.append(*args, **kwargs)
        return self.__class__(new_data)

    @property
    def ix(self):
        return self.data.ix

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
        if name in self.columns:
            return self[name]
        raise AttributeError("{0!r} has no attribute {1!r}".format(type(self).__name__, name))

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

        self_with_lens = self.data.join(self.lengths)

        new_regions_data = []

        for ix, row in self_with_lens.iterrows():
            remainder = row['length'] % resolution

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

            new_regions_data.append(row[['chromosome', 'start', 'end']])

        df = pd.DataFrame(new_regions_data, index=self.index)
        return Regions(df)

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
        ends = tss_locations + window_width

        regions_df = pd.DataFrame({'chromosome' : self['chromosome'],
                                   'start' : starts,
                                   'end' : ends}, index=self.index)

        return Regions(regions_df)

    @classmethod
    def from_encode_known_genes(cls, encode_known_genes_filename):
        from dgw.data.parsers import read_encode_known_genes
        return read_encode_known_genes(encode_known_genes_filename)

    @classmethod
    def from_gtf(cls, gtf_filename):
        from dgw.data.parsers import read_gtf
        return read_gtf(gtf_filename)

