__author__ = 'saulius'
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter

def plot_dtw_cost_and_path(cost_matrix, path, ax=None):
    '''
        Plots the DTW cost matrix and path.
        Code taken from the mlpy documentation:
        http://mlpy.sourceforge.net/docs/3.5/dtw.html
    :param cost_matrix:
    :param path:
    :return:
    '''

    if ax is None:
        ax = plt.gca()

    plot1 = ax.imshow(cost_matrix.T, origin='lower', cmap=cm.jet, aspect='auto')
    plt.colorbar(plot1)
    plot2 = ax.plot(path[0], path[1], 'w')
    xlim = ax.set_xlim((-0.5, cost_matrix.shape[0]-0.5))
    ylim = ax.set_ylim((-0.5, cost_matrix.shape[1]-0.5))

def plot_dtw_sequences_cost_and_path(sequence_x, sequence_y, cost_matrix, path):

    # Based on: http://matplotlib.org/examples/pylab_examples/scatter_hist.html

    null_fmt   = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02

    rect_cost = [left, bottom, width, height]

    rect_sequence_x = [left, bottom_h, width, 0.2]
    rect_sequence_y = [left_h, bottom, 0.2, height]

    ax_cost = plt.axes(rect_cost)
    ax_sequence_x = plt.axes(rect_sequence_x)
    ax_sequence_y = plt.axes(rect_sequence_y)

    # Drop labels
    ax_sequence_x.xaxis.set_major_formatter(null_fmt)
    ax_sequence_y.yaxis.set_major_formatter(null_fmt)

    # Plot main plot
    plot_dtw_cost_and_path(cost_matrix, path, ax = ax_cost)

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

    # rotate ax_sequence_y plot


    return ax_cost, ax_sequence_x, ax_sequence_y


def plot_dtw_path(path, ax = None, *args, **kwargs):
    if ax is None:
        ax = plt.gca()

    plot = ax.plot(path[0], path[1], *args, **kwargs)


