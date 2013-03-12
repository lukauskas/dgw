import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter

from distance import dtw_std

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

    # Based on: http://matplotlib.org/examples/pylab_examples/scatter_hist.html

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

    return ax_cost, ax_sequence_x, ax_sequence_y


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
    plot_dtw_sequences_dist_cost_and_path(sequence_x, sequence_y, dist, cost, path, *args, **kwargs)