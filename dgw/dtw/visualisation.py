from matplotlib.patches import ConnectionPatch
from dgw.util.plotting import pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter
import numpy as np
from dgw.dtw import reverse_sequence

from distance import dtw_std, dtw_path_is_reversed


def plot_dtw_cost_and_path(cost_matrix, path, ax=None):
    """
        Plots the DTW cost matrix and path.
        Code taken from the mlpy documentation:
        http://mlpy.sourceforge.net/docs/3.5/dtw.html
    :param cost_matrix:
    :param path:
    :return:
    """

    if ax is None:
        ax = plt.gca()

    plot1 = ax.imshow(cost_matrix.T, origin='lower', cmap=cm.gray, aspect='auto', interpolation='nearest')
    plt.colorbar(plot1)
    plot2 = ax.plot(path[0], path[1], 'w')

    xlim = ax.set_xlim((-0.5, cost_matrix.shape[0] - 0.5))
    ylim = ax.set_ylim((-0.5, cost_matrix.shape[1] - 0.5))

def plot_dtw_sequences_dist_cost_and_path(sequence_x, sequence_y, dist, cost_matrix, path):
    """
    Generates a full DTW plot with both the cost matrix, path and both of the sequences

    :param sequence_x:
    :param sequence_y:
    :param dist:
    :param cost_matrix:
    :param path:
    :return:
    """

    sequence_x = np.asarray(sequence_x)
    sequence_y = np.asarray(sequence_y)

    # Based on: http://matplotlib.org/examples/pylab_examples/scatter_hist.html
    figure = plt.figure()
    null_fmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    ax_width = 0.2
    bottom_h = left_h = left + width + 0.02

    rect_cost = [left, bottom, width, height]

    rect_sequence_x = [left, bottom_h, width, ax_width]
    rect_sequence_y = [left_h, bottom, ax_width, height]

    text_pos_x = left_h
    text_pos_y = bottom_h

    ax_cost = plt.axes(rect_cost)
    ax_sequence_x = plt.axes(rect_sequence_x)
    ax_sequence_y = plt.axes(rect_sequence_y)

    # Drop labels
    ax_sequence_x.xaxis.set_major_formatter(null_fmt)
    ax_sequence_y.yaxis.set_major_formatter(null_fmt)

    # Plot main plot
    plot_dtw_cost_and_path(cost_matrix, path, ax=ax_cost)

    ax_sequence_x.plot(range(len(sequence_x)), sequence_x)
    ax_sequence_y.plot(sequence_y, range(len(sequence_y)))

    # Set limits based on cost
    ax_sequence_x.set_xlim(ax_cost.get_xlim())
    ax_sequence_y.set_ylim(ax_cost.get_ylim())

    # Set Limits based on seq
    ax_sequence_x.set_ylim(min(sequence_x.min(), sequence_y.min()) - 1,
                           max(sequence_x.max(), sequence_y.max()) + 1)

    ax_sequence_y.set_xlim(min(sequence_x.min(), sequence_y.min()) - 1,
                           max(sequence_x.max(), sequence_y.max()) + 1)

    # Add test showing the value of dist
    plt.figtext(text_pos_x, text_pos_y, 'Distance:\n{0:.5f}'.format(dist), size='medium')

    return figure

def plot_dtw_path(path, ax=None, *args, **kwargs):
    """
    Plots path onto axis supplied by ax (or `plt.gca()` if no such axis given.

    :param path:
    :param ax:
    :param args:
    :param kwargs:
    :return:
    """
    if ax is None:
        ax = plt.gca()

    plot = ax.plot(path[0], path[1], *args, **kwargs)

def visualise_dtw(sequence_x, sequence_y, dtw_function=dtw_std, *args, **kwargs):
    """
    Visualises the DTW distance calculation between sequences x and y using the dtw fucntion provided
    
    :param sequence_x:
    :param sequence_y:
    :param dtw_function:
    :return:
    """
    dist, cost, path = dtw_function(sequence_x, sequence_y, dist_only=False)
    return plot_dtw_sequences_dist_cost_and_path(sequence_x, sequence_y, dist, cost, path, *args, **kwargs)

def visualise_dtw_mappings(sequence_x, sequence_y, dtw_function=dtw_std, columns=None, title=None, sequence_x_label=None,
                           sequence_y_label=None):

    def major_tick_step(ax, axis):
        if axis == 'x':
            ticks = ax.get_xticks()
        else:
            ticks = ax.get_yticks()

        try:
            step = ticks[1] - ticks[0]
        except IndexError:
            step = 0.2

        return step

    def expand_axes(ax):
        x_increment = major_tick_step(ax, 'x') / 8.0

        min_x, max_x = ax.get_xlim()
        ax.set_xlim(min_x - x_increment, max_x + x_increment)

        y_increment = major_tick_step(ax, 'y') / 8.0
        min_y, max_y = ax.get_ylim()
        ax.set_ylim(min_y - y_increment, max_y + y_increment)

    def add_reversed_annotation(ax):
        min_x, max_x = ax.get_xlim()
        min_y, max_y = ax.get_ylim()

        offset_x = major_tick_step(ax, 'x')
        offset_y = major_tick_step(ax, 'y')
        ax.text(min_x + offset_x / 8, max_y - offset_y / 2 - offset_y / 8, '(reversed)')

    dist, cost, path = dtw_function(sequence_x, sequence_y, dist_only=False)

    sequence_x = np.asarray(sequence_x)
    sequence_y = np.asarray(sequence_y)

    try:
        ndim = sequence_x.shape[1]
    except IndexError:
        ndim = 1

    reversed = dtw_path_is_reversed(path)
    if reversed:
        sequence_x = reverse_sequence(sequence_x)
        path_x = np.max(path[0]) - path[0]
        path_y = path[1]
        path = (path_x, path_y)

    sequence_y_T = np.atleast_2d(sequence_y.T)
    sequence_x_T = np.atleast_2d(sequence_x.T)

    if columns is None and ndim > 1:
        columns = ['Dimension #{0}'.format(i) for i in range(1, ndim+1)]
    elif ndim > 1:
        if len(columns) != ndim:
            raise ValueError('Number of column titles does not match the number of dimensions')

    main_y_axis = None
    xaxes_regular = [None] * ndim
    xaxes_warped = [None] * ndim
    figure = plt.figure()
    figure.subplots_adjust(wspace=0.01, hspace=0.1)
    for i in range(ndim):
        x = sequence_x_T[i]
        print x.shape
        y = sequence_y_T[i]


        ax2 = plt.subplot(2, ndim, ndim + i + 1, sharey=main_y_axis, sharex=xaxes_warped[i])
        ax2.plot(y, color='g')

        if i > 0:
            ax2.yaxis.set_visible(False)

        expand_axes(ax2)

        if not main_y_axis:
            main_y_axis = ax2

        if not xaxes_warped[i]:
            xaxes_warped[i] = ax2

        ax1 = plt.subplot(2, ndim, i + 1, sharey=main_y_axis, sharex=xaxes_regular[i])
        ax1.plot(x, color='b')
        expand_axes(ax1)

        if ndim > 1:
            ax1.set_title(columns[i])

        if i > 0:
            ax1.yaxis.set_visible(False)


        if not xaxes_regular[i]:
            xaxes_regular[i] = ax1
        if reversed:
            add_reversed_annotation(ax1)

        for p_i, p_j in zip(path[0], path[1]):
            xy_a = (p_i, x[p_i])
            xy_b = (p_j, y[p_j])
            con = ConnectionPatch(xyA=xy_a, xyB=xy_b, coordsA="data", coordsB="data",
                                 axesA=ax1, axesB=ax2, arrowstyle="-", shrinkB=2, shrinkA=2, alpha=0.2)

            ax1.add_artist(con)

    if title is not None:
        plt.suptitle(title)

    lines = xaxes_regular[0].get_lines()
    lines.extend(xaxes_warped[0].get_lines())

    if sequence_x_label and sequence_y_label:
        plt.figlegend(lines, (sequence_x_label, sequence_y_label), 'lower center')

    return figure