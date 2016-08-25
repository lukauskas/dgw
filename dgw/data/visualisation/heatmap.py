from collections import defaultdict
from logging import debug
from math import ceil
from dgw.util.plotting import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FuncFormatter, IndexLocator
import palettable
import numpy as np
from dgw.util.plotting import matplotlib
import pandas as pd
from dgw.dtw.utilities import no_nans_len


def dataset_ticks(dataset, scale=1):
    """
    Ticks formater for heatmap. Shows dataset indices as ticks.

    Use it as `ax.yaxis.set_major_formatter(dataset_ticks(dataset, scale))`

    :param dataset:
    :param scale:
    :return:
    """
    if scale is None:
        scale = 1

    def f(tick, pos):
        try:
            return dataset.items[int((tick / scale) + 0.5)]
        except IndexError:
            return ''

    return FuncFormatter(f)

def raw_plot_data_as_heatmap(data_frame, ax=None, highlight_masks=None, highlight_colours=None, *args, **kwargs):
    if ax is None:
        ax = plt.gca()

    # Create a masked array from data_frame
    masked_arr = np.ma.array(data_frame, mask=np.isnan(data_frame))

    # Use colormap from colorbrewer, draw NaNs as white
    cmap = palettable.colorbrewer.sequential.YlGnBu_9.mpl_colormap

    cmap.set_bad('#FFFFFF', 1.0)
    HEATMAP_COLORS = ['k', 'w']
    result = ax.imshow(masked_arr, origin='lower', cmap=cmap, aspect="auto", interpolation='nearest', *args, **kwargs)
    if highlight_masks is not None:
        for j, mask in highlight_masks.iteritems():
            highlight_cmap = matplotlib.colors.ListedColormap([highlight_colours[j]])
            highlight_cmap.set_bad('#000000', 0)
            ax.imshow(highlight_masks[j], origin='lower', cmap=highlight_cmap, aspect="auto", interpolation='nearest',
                               *args, **kwargs)

    return result

def plot(alignments, titles=None, horizontal_grid=True,
         no_y_axis=False, sort_by=None, subplot_spec=None, share_y_axis=None, scale_y_axis=None, highlighted_points={},
         highlight_colours=None,
         rasterized=True, replace_with_dataframe=None):
    """

    :param alignments: `AlignmentsData` object
    :param titles: Titles of the heatmaps. If set to none `alignments.dataset_axis` will be used
    :param horizontal_grid: Whether to plot heatmap on a horizontal grid
    :param no_y_axis: Will not plot the major (items) axis.
    :param sort_by: Sort values by 'length' or supplied index.
    :param subplot_spec: SubplotSpec of the suplot to plot hte heatmaps in. See `matplotlib.gridspec` package.
    :param share_y_axis: if not None, the plot will share the major (items) axis with the specified axis
    :param scale_y_axis: plot will scale the y axis by the specified number (linearly) if set.
    :type scale_y_axis: int
    :param highlighted_points: a (index, array) dictionary of points that should be highlighted in the heatmap
    :param: rasterized: whether to rasterize the plot or not (faster rending for rasterized)
    :param replace_with_dataframe: replace the data with the dataframe provided. Use for toggleable overlays
    :return: returns the shared y axis
    """
    from dgw.util.plotting import pyplot as plt
    from matplotlib import gridspec

    number_of_datasets = alignments.number_of_datasets

    if titles is None:
        titles = alignments.dataset_axis

    # Should be OK just to just sort the first dataset as all others should have the same len
    sample_data_frame = alignments.dataset_xs(alignments.dataset_axis[0], copy=False).T
    lengths = sample_data_frame.apply(no_nans_len, axis=1)
    max_len = lengths.max()

    # Sorting
    if isinstance(sort_by, pd.Index):
        sorted_index = sort_by
    elif sort_by == 'length':
        debug('Sorting by length')
        lengths.sort()  # Should do in-place
        sorted_index = lengths.index
    elif sort_by is None:
        sorted_index = alignments.items
    else:
        raise ValueError('Unsupported sort_by value provided: {0!r}. Only None or \'length\' supported'.format(sort_by))
    # Apply sorting
    alignments = alignments.ix[sorted_index]
    # Create the instance of axis formatter
    tick_formatter = dataset_ticks(alignments, scale_y_axis)

    # Subplot creation

    if horizontal_grid:
        grid = (1, number_of_datasets + 1)
        width_ratios = [5] * number_of_datasets
        width_ratios.append(1)
        spacing_kwargs = {'wspace': 0.01, 'width_ratios': width_ratios}  # Almost no space between plots
    else:
        grid = (number_of_datasets, 2)
        spacing_kwargs = {'hspace': 0.15, 'width_ratios': [5, 1]} # Allow a bit of a space for title

    if not subplot_spec:
        gs = gridspec.GridSpec(*grid, **spacing_kwargs)
    else:
        gs = gridspec.GridSpecFromSubplotSpec(*grid, subplot_spec=subplot_spec, **spacing_kwargs)

    # Main drawing loop
    first_axis = None
    extent = None
    if scale_y_axis:
        extent = [0, max_len, 0, alignments.number_of_items * scale_y_axis] # Multiply by 10 as that is what matplotlib's dendrogram returns

    if highlighted_points:
        highlight_masks = defaultdict(lambda: np.zeros((alignments.number_of_items, max_len), dtype=np.bool))
        for i, ix in enumerate(sorted_index):
            try:
                points_of_interest = highlighted_points[ix]
            except KeyError:
                continue

            for j, points in points_of_interest.iteritems():
                highlight_masks[j][i][points] = 1

        for j in highlight_masks.iterkeys():
            highlight_masks[j] = np.ma.masked_where(highlight_masks[j] <= 0, highlight_masks[j])

    else:
        highlight_masks = None

    for i, (ix, title) in enumerate(zip(alignments.dataset_axis, titles)):
        t_gs = gs[:, i] if horizontal_grid else gs[i, 1]
        if i == 0:
            if not share_y_axis:
                first_axis = plt.subplot(t_gs)
                share_y_axis = first_axis
            else:
                first_axis = plt.subplot(t_gs, sharey=share_y_axis)
        else:
            # Remember to share axes
            plt.subplot(t_gs, sharex=first_axis, sharey=share_y_axis)

        plt.gca().get_yaxis().set_major_formatter(tick_formatter)
        #plt.gca().get_yaxis().set_major_locator(IndexLocator(1000, 0))
        # Remove redundant axes
        if horizontal_grid:
            if i > 0 or no_y_axis:
                plt.gca().get_yaxis().set_visible(False)
        else:
            # Leave only last axis
            if i + 1 < number_of_datasets:
                plt.gca().get_xaxis().set_visible(False)

            if no_y_axis:
                plt.gca().get_yaxis().set_visible(False)

        data_to_plot = alignments.dataset_xs(ix, copy=False).T
        # Cut most of the columns that are NaNs out of the plot
        data_to_plot = data_to_plot[data_to_plot.columns[:max_len]]

        if replace_with_dataframe is not None:
            result = raw_plot_data_as_heatmap(replace_with_dataframe, extent=extent,
                                              rasterized=rasterized)
        else:
            result = raw_plot_data_as_heatmap(data_to_plot, extent=extent,
                                          highlight_masks=highlight_masks, rasterized=rasterized, highlight_colours=highlight_colours)
        debug(plt.gca().get_ylim())
        debug(plt.gca().get_xlim())

        plt.title(title)
        plt.gca().title.set_fontsize(7)
        plt.gca().grid(False)
        plt.setp(plt.gca().get_xticklabels(), rotation='vertical', fontsize=7)

    # Colorbar
    colorbar_axis = plt.subplot(gs[:, number_of_datasets] if horizontal_grid else gs[:, 2])
    plt.colorbar(result, cax=colorbar_axis)

    return first_axis