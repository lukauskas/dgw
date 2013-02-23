import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import visualisation.heatmap as heatmap

class AlignmentsDataIndexer(object):
    """
    A wrapper around `_NDFrameIndexer` that would return `AlignmentsData` objects instead of `pd.Panel` objects
    """
    _ndframe_indexer = None

    def __init__(self, ndframe_indexer):
        self._ndframe_indexer = ndframe_indexer

    def __getitem__(self, key):
        result = self._ndframe_indexer.__getitem__(key)
        if isinstance(result, pd.Panel):
            return AlignmentsData(result)
        else:
            return result

    def __setitem__(self, key, value):
        self._ndframe_indexer.__setitem__(self, key, value)

class AlignmentsData(object):

    _data = None
    _scale = None
    def __init__(self, panel, scale='raw'):
        """
        Initialises `AlignmentsData` with a `panel` provided.
        The panel is assumed to have data sets on the minor axis
        See `dgw.data.parsers.read_bam` for how to generate this data.

        :param panel: `pd.Panel` object `AlignmentsData` will be initialised on, or `pd.DataFrame` that will be converted
                     to Panel
        :param scale: the scale of data
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

        self._scale = scale

    @property
    def data(self):
        return self._data

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
    def dataset_axis(self):
        return self.data.minor_axis

    def dataset_xs(self, *args, **kwargs):
        return self.data.minor_xs(*args, **kwargs)

    @property
    def items(self):
        return self.data.items

    @property
    def ix(self):
        return AlignmentsDataIndexer(self.data.ix)

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

        new_data = (self.data + 1).apply(np.log)
        return AlignmentsData(new_data, scale='log')

    def plot_heatmap(self, titles=None, horizontal_grid=True, parent_subplot_spec=None, *args, **kwargs):
        """
        Plots heatmap of the data stored in the panel.

        :param args: args to be passed into `data.visualisation.heatmap.plot`
        :param kwargs: kwargs to be passed into `data.visualisation.heatmap.plot`
        :return:
        """
        import matplotlib.pyplot as plt
        from matplotlib import gridspec

        number_of_datasets = self.number_of_datasets
        if titles is None:
            titles = self.dataset_axis

        if number_of_datasets > 1:
            if horizontal_grid:
                grid = (1, number_of_datasets)
            else:
                grid = (number_of_datasets, 1)

            if not parent_subplot_spec:
                gs = gridspec.GridSpec(*grid)
            else:
                gs = gridspec.GridSpecFromSubplotSpec(*grid, subplot_spec=parent_subplot_spec)

        for i, (ix, title) in enumerate(zip(self.dataset_axis, titles)):
            if number_of_datasets > 1:
                t_gs = gs[i]
                plt.subplot(t_gs)
            elif parent_subplot_spec:
                plt.subplot(parent_subplot_spec)

            data_to_plot = self.dataset_xs(ix, copy=False).T
            heatmap.plot(data_to_plot, *args, **kwargs)
            plt.title(title)

    def plot_all(self):
        for ix in self.items:
            plt.plot(self.ix[ix], label=ix)


    def plot_random_sample(self, k):
        import matplotlib.pyplot as plt
        import random

        items = self.items
        for i in xrange(k):
            ax = plt.subplot(1, k, i+1)
            random_item = self.ix[random.choice(items)]
            random_item.plot(ax=ax)




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

    def to_bed(self, bed_file):
        from dgw.data.parsers import write_bed
        return write_bed(self, bed_file)

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

    def regions_not_in_dataset(self, dataset):
        """
        Returns all regions that are not in the dataset provided.

        :param dataset:
        :type dataset: AlignmentsData
        :return:
        """
        missing_indices = self.index[~self.index.isin(dataset.items)]
        return Regions(self.data.ix[missing_indices])



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

