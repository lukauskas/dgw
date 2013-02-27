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
        # A5 Paper size
        plt.figure(num=None, figsize=(8.3, 5.8), facecolor='w', edgecolor='k')
        self.title = plt.suptitle("") # Create text object for title
        self.draw_buttons()
        self.draw()
        plt.show()

class HierarchicalClusteringViewer(object):
    _hierarchical_clustering_object = None
    _gs_main = None

    def __init__(self, hierarchical_clustering_object):
        self._hierarchical_clustering_object = hierarchical_clustering_object

    @property
    def gs_dendrogram(self):
        return self.gs_main[0]

    @property
    def gs_heatmap(self):
        return self.gs_main[1]

    @property
    def gs_main(self):
        return self._gs_main

    @property
    def hierarchical_clustering_object(self):
        return self._hierarchical_clustering_object

    def set_up(self):
        self._gs_main = gridspec.GridSpec(1, 2, wspace=0)

    def draw(self):
        ax_dendrogram = plt.subplot(self.gs_dendrogram)
        hc = self.hierarchical_clustering_object

        dendrogram_dict = hc.dendrogram(orientation='right', get_leaves=True,
                                        color_threshold=0, distance_sort=True)
        leaves = dendrogram_dict['leaves']

        index = hc.data.items[leaves]

        DENDROGRAM_SCALE = 10  # scipy.cluster.hierarachy.dendrogram scales all y axis values by tenfold for some reason
        hc.data.plot_heatmap(subplot_spec=self.gs_heatmap, no_y_axis=True, sort_by=index, share_y_axis=ax_dendrogram,
                             scale_y_axis=DENDROGRAM_SCALE)

    def show(self):
        self.set_up()
        self.draw()
        plt.show()
