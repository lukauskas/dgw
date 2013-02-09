__author__ = 'saulius'
import pandas as pd
import data.visualisation.heatmap as heatmap
import matplotlib.pyplot as plt

class AggregatedAlignmentsPanel(pd.Panel):
    def mean(self, axis='items', skipna=True):
        # Override axis parameter in the pd.Panel mean function
        return super(AggregatedAlignmentsPanel, self).mean(axis=axis, skipna=skipna)

    @property
    def number_of_datasets(self):
        return len(self.minor_axis)

    def plot_heatmap(self, *args, **kwargs):
        """
        Plots heatmap of the data stored in the panel.

        :param args: args to be passed into `data.visualisation.heatmap.plot`
        :param kwargs: kwargs to be passed into `data.visualisation.heatmap.plot`
        :return:
        """
        number_of_datasets = self.number_of_datasets

        for i, title in enumerate(self.minor_axis):
            if number_of_datasets > 1:
                plt.subplot(1, number_of_datasets, i+1) # TODO: consider doing sublot with multiple lines

            data_to_plot = self.minor_xs(title, copy=False).T
            heatmap.plot(data_to_plot, *args, **kwargs)
            plt.title(title)

class Regions(pd.DataFrame):
    REQUIRED_COLUMNS = frozenset(['chromosome', 'start', 'end'])

    def __init__(self, *args, **kwargs):
        super(Regions, self).__init__(*args, **kwargs)

        # Verify that all required columns are in the DF
        for column in self.REQUIRED_COLUMNS:
            if column not in self.columns:
                raise ValueError('No such column {0!r} in provided DataFrame'.format(column))

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
        Adjusts the region boundaries to fit the resolution by extending the regions boundaries as to length of the regions
        is divisible from resolution.

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

        self_with_lens = self.join(self.lengths)

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

            row['end']   += add_right

            new_regions_data.append(row[['chromosome', 'start', 'end']])

        new_regions = Regions(new_regions_data, index=self.index)
        return new_regions
