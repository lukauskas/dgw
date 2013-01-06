__author__ = 'saulius'
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def plot_dtw_cost_and_path(cost_matrix, path):
    '''
        Plots the DTW cost matrix and path.
        Code taken from the mlpy documentation:
        http://mlpy.sourceforge.net/docs/3.5/dtw.html
    :param cost_matrix:
    :param path:
    :return:
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot1 = plt.imshow(cost_matrix.T, origin='lower', cmap=cm.gray, interpolation='nearest')
    plot2 = plt.plot(path[0], path[1], 'w')
    xlim = ax.set_xlim((-0.5, cost_matrix.shape[0]-0.5))
    ylim = ax.set_ylim((-0.5, cost_matrix.shape[1]-0.5))


def plot_dtw_path(path, *args, **kwargs):

    plot = plt.plot(path[0], path[1], *args, **kwargs)


