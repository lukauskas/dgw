from logging import debug
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button


class ClusterPreviewer(object):
    _clusters = None
    _main_gs = None
    _current_cluster_id = None

    def __init__(self, cluster_roots):
        self._clusters = cluster_roots
        self._initialise_grid()

        self._current_cluster_id = 0

        debug('Calculating projections')
        for c in self._clusters:
            c.ensure_projections_are_calculated()


    def gs_prototype(self):
        return self._main_gs[1, 0]
    def gs_projected_mean(self):
        return self._main_gs[1, 1]
    def gs_projected_heatmap(self):
        return self._main_gs[3, :]
    def gs_heatmap(self):
        return self._main_gs[2, :]


    def _initialise_grid(self):
        self._main_gs = gridspec.GridSpec(5, 2, height_ratios=[1, 40, 20, 20, 4])


    def current_cluster(self):
        return self._clusters[self._current_cluster_id]

    def _callback_next(self, event):
        debug('next clicked')
        self._current_cluster_id = (self._current_cluster_id + 1) % len(self._clusters)
        self.draw()

    def _callback_previous(self, event):
        debug('previous clicked')
        self._current_cluster_id = (self._current_cluster_id - 1) % len(self._clusters)
        self.draw()

    def draw_buttons(self):
        buttons_left = 0.05
        buttons_bottom = 1 - 0.1

        button_width = 0.1
        button_height = 0.05

        button_spacing = 0.02

        ax_button_prev = plt.axes([buttons_left, buttons_bottom, button_width, button_height])
        self.button_prev = Button(ax_button_prev, 'Previous')
        self.button_prev.on_clicked(self._callback_previous)

        ax_button_next = plt.axes([buttons_left + button_width + button_spacing, buttons_bottom,
                                   button_width, button_height])
        self.button_next = Button(ax_button_next, 'Next')
        self.button_next.on_clicked(self._callback_next)

    def draw(self):

        ax_prototype = plt.subplot(self.gs_prototype())
        plt.cla()
        plt.title('Prototype')
        current_cluster = self.current_cluster()
        current_cluster.prototype.plot(ax=ax_prototype, legend=False)

        plt.figlegend(*ax_prototype.get_legend_handles_labels(), loc='lower center')

        self.title.set_text('Cluster #{0}/{2} ({1} elements)'.format(self._current_cluster_id + 1,
                                                                     current_cluster.n_items,
                                                                     len(self._clusters)))

        plt.subplot(self.gs_heatmap())
        plt.cla()
        current_cluster.data.plot_heatmap(horizontal_grid=True, subplot_spec=self.gs_heatmap())

        # Projections
        projections = current_cluster.projected_data

        ax_projections = plt.subplot(self.gs_projected_mean())
        plt.cla()
        plt.title('Average DTW projection onto prototype')
        projections.mean().plot(ax=ax_projections, legend=False)

        ax_projected_heatmap = plt.subplot(self.gs_projected_heatmap())
        plt.cla()
        current_cluster.projected_data.plot_heatmap(horizontal_grid=True, subplot_spec=self.gs_projected_heatmap())

        # Finally issue a draw command for the plot
        plt.draw()


    def show(self):
        # A4 Paper size
        plt.figure(num=None, figsize=(11.7, 8.3), dpi=100, facecolor='w', edgecolor='k')
        self.title = plt.suptitle("") # Create text object for title
        self.draw_buttons()
        self.draw()
        plt.show()

class HierarchicalClusteringViewer(object):
    _hierarchical_clustering_object = None
    _gs_main = None
    _figure = None
    _ax_dendrogram = None
    _line = None
    _cut_xdata = 0

    def __init__(self, hierarchical_clustering_object):
        self._hierarchical_clustering_object = hierarchical_clustering_object

    @property
    def gs_dendrogram(self):
        return self.gs_main[0]

    @property
    def gs_heatmap(self):
        return self.gs_main[1]

    @property
    def ax_dendrogram(self):
        return self._ax_dendrogram

    @property
    def gs_main(self):
        return self._gs_main

    @property
    def hierarchical_clustering_object(self):
        return self._hierarchical_clustering_object

    def set_up(self):
        self._gs_main = gridspec.GridSpec(1, 2, wspace=0)
        self._figure = plt.gcf()
        self._ax_dendrogram = plt.subplot(self.gs_dendrogram)
        self._figure.canvas.mpl_connect('button_press_event', self._onclick_listener)

    def draw_dendrogram(self):
        ax_dendrogram = self.ax_dendrogram
        plt.subplot(ax_dendrogram)
        hc = self.hierarchical_clustering_object

        dendrogram_dict = hc.dendrogram(orientation='right', get_leaves=True,
                                        distance_sort=True, color_threshold=self._cut_xdata)
        leaves = dendrogram_dict['leaves']

        return leaves

    def draw(self):
        leaves = self.draw_dendrogram()
        hc = self.hierarchical_clustering_object
        index = hc.data.items[leaves]

        DENDROGRAM_SCALE = 10  # scipy.cluster.hierarachy.dendrogram scales all y axis values by tenfold for some reason
        hc.data.plot_heatmap(subplot_spec=self.gs_heatmap, no_y_axis=True, sort_by=index, share_y_axis=self.ax_dendrogram,
                             scale_y_axis=DENDROGRAM_SCALE)

    def show(self):
        # A5 Paper size
        plt.figure(num=None, figsize=(8.3, 5.8), dpi=100, facecolor='w', edgecolor='k')
        self.set_up()
        self.draw()
        plt.show()

    def cut_line(self, xdata):
        if self._line is not None:
            self._line.set_xdata(xdata)
        else:
            self._line = self._ax_dendrogram.axvline(linestyle='--', x=xdata, color='k')

        self._cut_xdata = xdata
        self.draw_dendrogram()


    def _onclick_listener(self, event):
        # Allow only in axis, only the left button and only regions with value
        if event.inaxes != self.ax_dendrogram or event.button != 1 or event.xdata is None:
            return
        self.cut_line(event.xdata)

