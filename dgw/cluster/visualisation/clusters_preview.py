import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button


class ClusterAssignmentsPreviewer(object):
    _clusters = None
    _main_gs = None
    _current_cluster_id = None

    def __init__(self, cluster_assignments):
        self._clusters = cluster_assignments
        self._initialise_grid()

        self._current_cluster_id = 0

    def gs_prototype(self):
        return plt.subplot(self._main_gs[1, 0])
    def ax_prototype(self):
        return plt.subplot(self.gs_prototype())

    def gs_heatmap(self):
        return self._main_gs[1, 1]
    def ax_heatmap(self):
        return plt.subplot(self.gs_heatmap())

    def _initialise_grid(self):
        self._main_gs = gridspec.GridSpec(3, 2, height_ratios=[1, 40, 4])


    def current_cluster(self):
        return self._clusters[self._current_cluster_id]

    def draw_buttons(self):
        buttons_left = 0.05
        buttons_top = 1 - 0.05

        button_width = 0.08
        button_height = 0.04

        button_spacing = 0.02

        ax_button_prev = plt.axes([buttons_left, buttons_top, button_width, button_height])
        button_prev = Button(ax_button_prev, 'Previous')

        ax_button_next = plt.axes([buttons_left, buttons_top - button_height - button_spacing,
                                   button_width, button_height])
        button_next = Button(ax_button_next, 'Next')

    def draw(self):

        ax_prototype = self.ax_prototype()
        current_cluster = self.current_cluster()
        current_cluster.prototype.plot(ax=ax_prototype, legend=False)

        plt.figlegend(*ax_prototype.get_legend_handles_labels(), loc='lower center')

        plt.suptitle('Cluster #{0} ({1} elements)'.format(self._current_cluster_id, current_cluster.n_items))

        ax_heatmap = self.ax_heatmap()
        current_cluster.items.plot_heatmap(horizontal_grid=False, parent_subplot_spec=self.gs_heatmap())

        self.draw_buttons()


    def show(self):
        # A5 Paper size
        plt.figure(num=None, figsize=(8.3, 5.8), dpi=150, facecolor='w', edgecolor='k')
        self.draw()
        plt.show()
