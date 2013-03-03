from logging import debug
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button
import Tkinter, tkFileDialog

class ClusterPreviewer(object):
    _clusters = None
    _main_gs = None
    _current_cluster_id = None

    root_window = None
    _add_save_button = None
    def __init__(self, cluster_roots):
        self._clusters = cluster_roots
        self._initialise_grid()

        self._current_cluster_id = 0

        debug('Calculating projections')
        self._add_save_button = True
        for c in self._clusters:
            if c.regions is None:
                debug('Cluster {0!r} has no regions information, cannot save'.format(c))
                self._add_save_button = False

            c.ensure_projections_are_calculated()
            c.ensure_points_of_interest_are_tracked_down()


    def gs_prototype(self):
        return self._main_gs[1, 0]
    def gs_projected_mean(self):
        return self._main_gs[1, 1]
    def gs_projected_heatmap(self):
        return self._main_gs[3, :]
    def gs_heatmap(self):
        return self._main_gs[2, :]


    def _initialise_grid(self):
        self._main_gs = gridspec.GridSpec(5, 2, height_ratios=[1, 40, 20, 20, 4], hspace=0.4)


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

    def _callback_save(self, event):
        debug('Save clicked')

        # Create a root window that can be closed as otherwise tkFileDialog will create unclosable one
        if self.root_window is None:
            self.root_window = Tkinter.Tk()
            self.root_window.withdraw()  # Closing the newly-created win

        directory = tkFileDialog.askdirectory(title="Select a directory where the resulting clusters will be saved to")
        if not directory: # Filedialog returns empty string if escape is pressed
            return

        if not os.path.exists(directory):
            os.makedirs(directory)

        print '> Saving clusters to directory {0}'.format(directory)
        for i, c in enumerate(self._clusters):
            filename = os.path.join(directory, 'dgw_cluster_{0}.bed'.format(i + 1))
            c.save_as_encode_track(filename, track_name='DGWCluster{0}'.format(i + 1),
                                   track_description='DGW Cluster {0}/{1}'.format(i + 1, len(self._clusters)))
        print '> Saved'

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

        if self._add_save_button:
            ax_button_save = plt.axes([buttons_left + button_width * 2 + button_spacing * 2, buttons_bottom,
                                       button_width, button_height])
            self.button_save = Button(ax_button_save, 'Save all')
            self.button_save.on_clicked(self._callback_save)

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

        shared_axis = current_cluster.data.plot_heatmap(horizontal_grid=True, subplot_spec=self.gs_heatmap(),
                                                        sort_by=None, highlighted_points=current_cluster.points_of_interest)

        # Projections
        projections = current_cluster.projected_data

        ax_projections = plt.subplot(self.gs_projected_mean())
        plt.cla()
        plt.title('Average DTW projection onto prototype')
        projections.mean().plot(ax=ax_projections, legend=False)

        ax_projected_heatmap = plt.subplot(self.gs_projected_heatmap())
        plt.cla()
        current_cluster.projected_data.plot_heatmap(horizontal_grid=True, subplot_spec=self.gs_projected_heatmap(),
                                                    share_y_axis=shared_axis, sort_by=None, highlighted_points=current_cluster.tracked_points_of_interest)

        # Finally issue a draw command for the plot
        plt.draw()


    def show(self):
        # A5 Paper size
        plt.figure(num=None, figsize=(12, 10), facecolor='w', edgecolor='k')
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
    def clusters(self):
        if self._cut_xdata > 0:
            return self._hierarchical_clustering_object.cut(self._cut_xdata)
        else:
            return []

    @property
    def gs_dendrogram(self):
        return self.gs_main[1, 0]

    @property
    def gs_heatmap(self):
        return self.gs_main[1, 1]

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
        self._gs_main = gridspec.GridSpec(2, 2, wspace=0, height_ratios=[1, 15])
        self._figure = plt.gcf()
        self._ax_dendrogram = plt.subplot(self.gs_dendrogram)
        self._figure.canvas.mpl_connect('button_press_event', self._onclick_listener)

        self.draw_buttons()

    def draw_dendrogram(self):
        ax_dendrogram = self.ax_dendrogram
        plt.subplot(ax_dendrogram)
        hc = self.hierarchical_clustering_object

        dendrogram_dict = hc.dendrogram(orientation='right', get_leaves=True,
                                        distance_sort=True, color_threshold=self._cut_xdata)
        leaves = dendrogram_dict['leaves']
        plt.setp(plt.gca().get_xticklabels(), rotation='vertical', fontsize=7)
        plt.gca().get_yaxis().set_visible(False)
        return leaves

    def _callback_preview(self, event):
        if not self.clusters:
            return

        # Else, if we have clusters
        pw = ClusterPreviewer(self.clusters)
        pw.show()

    def _callback_save(self, event):
        pass

    def draw_buttons(self):
        buttons_left = 0.05
        buttons_bottom = 1 - 0.1

        button_width = 0.1
        button_height = 0.05

        button_spacing = 0.02

        ax_button_preview = plt.axes([buttons_left, buttons_bottom, button_width, button_height])
        self.button_preview = Button(ax_button_preview, 'Preview')
        self.button_preview.on_clicked(self._callback_preview)


    def draw(self):
        leaves = self.draw_dendrogram()
        hc = self.hierarchical_clustering_object
        index = hc.data.items[leaves]

        DENDROGRAM_SCALE = 10  # scipy.cluster.hierarachy.dendrogram scales all y axis values by tenfold for some reason
        hc.data.plot_heatmap(subplot_spec=self.gs_heatmap, no_y_axis=True, sort_by=index, share_y_axis=self.ax_dendrogram,
                             scale_y_axis=DENDROGRAM_SCALE, highlighted_points=hc.data.points_of_interest)

    def show(self):
        # A5 Paper size
        plt.figure(num=None, figsize=(12, 10), facecolor='w', edgecolor='k')
        self.set_up()
        self.draw()
        plt.show()

    def cut_line(self, xdata):

        self._cut_xdata = xdata
        self.draw_dendrogram()
        if self._line is not None:
            self._line.set_xdata(xdata)
        else:
            self._line = self._ax_dendrogram.axvline(linestyle='--', x=xdata, color='k')

        self._figure.canvas.draw()


    def _onclick_listener(self, event):
        # Allow only in axis, only the left button and only regions with value
        if event.inaxes != self.ax_dendrogram or event.button != 1 or event.xdata is None:
            return
        self.cut_line(event.xdata)

