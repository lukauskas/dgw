import pandas as pd
import matplotlib.pyplot as plt

import dgw.data.visualisation.heatmap as heatmap
from dgw.data.parsers import read_bed, read_bam

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
        return AlignmentsData(read_bam(bam_files))

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

    # --- Initialisation ----------------------------------------------------------------------------------------------
    @classmethod
    def from_bed(cls, bed_file):
        return cls(read_bed(bed_file))

    # --- Functions that provide direct access to the DataFrame behind all this ----------------------------------------
    def __getitem__(self, item):
        return self.data[item]

    def iterrows(self):
        return self.data.iterrows()

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


