from collections import defaultdict
from logging import debug, error
import os

from matplotlib.patches import Rectangle
from dgw.util.plotting import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button
import numpy as np
import pandas as pd

from dgw.dtw.visualisation import visualise_dtw_mappings


class ClusterPreviewer(object):
    _clusters = None
    _main_gs = None
    _current_cluster_id = None

    root_window = None
    _buttons = None

    _ax_heatmap = None
    _ax_projected_heatmap = None
    _enlarged_heatmap_axis = None
    _enlarged_projected_heatmap_axis = None


    def __init__(self, cluster_roots, output_directory, configuration_filename, cut_level, dtw_function, highlight_colours=None):
        self.output_directory = output_directory
        self.cut_level = cut_level
        self.configuration_filename = configuration_filename

        self._clusters = cluster_roots
        self._initialise_grid()

        self._current_cluster_id = 0

        self.dtw_function = dtw_function

        debug('Calculating projections')
        for c in self._clusters:

            c.ensure_projections_are_calculated()
            c.ensure_points_of_interest_are_tracked_down()

        self._buttons = []
        self.highlight_colours = highlight_colours

        self._warping_conservation_view = False

    def gs_prototype(self):
        return self._main_gs[1, 0]
    def gs_projected_mean(self):
        return self._main_gs[1:3, 1]
    def gs_projected_heatmap(self):
        return self._main_gs[4, :]
    def gs_heatmap(self):
        return self._main_gs[3, :]

    def gs_warping_preservation(self):
        return self._main_gs[2, 0]


    def ax_heatmap(self):
        if self._ax_heatmap is None:
            self._ax_heatmap = plt.subplot(self.gs_heatmap())
        return self._ax_heatmap

    def ax_projected_heatmap(self):
        if self._ax_projected_heatmap is None:
            self._ax_projected_heatmap = plt.subplot(self.gs_projected_heatmap())
        return self._ax_projected_heatmap

    def _initialise_grid(self):
        self._main_gs = gridspec.GridSpec(6, 2, height_ratios=[1, 30, 10, 20, 20, 4], hspace=0.4)


    def current_cluster(self):
        """

        :return:
        :rtype: DTWClusterNode
        """
        return self._clusters[self._current_cluster_id]

    def change_cluster(self, new_cluster_id):
        self._current_cluster_id = new_cluster_id
        self._warping_conservation_view = False


    def _callback_next(self, event):
        debug('next clicked')
        self.change_cluster((self._current_cluster_id + 1) % len(self._clusters))
        self.draw()

    def _callback_previous(self, event):
        debug('previous clicked')
        self.change_cluster((self._current_cluster_id - 1) % len(self._clusters))
        self.draw()

    def save_previewer_windows(self, output_directory):

        if not os.path.isdir(output_directory):
            os.makedirs(output_directory)

        self.create_figure(interactive=False)

        for i in range(len(self._clusters)):
            if i > 0:
                self.change_cluster((self._current_cluster_id + 1) % len(self._clusters))

            cluster_id = self._current_cluster_id
            filename = 'summary-{}.pdf'.format(cluster_id+1)

            self.draw()
            plt.savefig(os.path.join(output_directory, filename),
                       create=False)


    def save_clusters(self, directory):

        directory = self.output_directory
        if not os.path.exists(directory):
            os.makedirs(directory)

        repro_filename = os.path.join(directory, 'repro.txt')
        repro_file = open(repro_filename, 'w')

        try:
            repro_file.write('To reproduce the cut, '
                             'pass the following parameter to dgw-explorer:\n-c {0}\n'.format(self.cut_level))
        finally:
            repro_file.close()

        bed_directory = os.path.join(directory, 'items')
        if not os.path.exists(bed_directory):
            os.makedirs(bed_directory)

        heatmaps_directory = os.path.join(directory, 'heatmaps')
        if not os.path.exists(heatmaps_directory):
            os.makedirs(heatmaps_directory)

        prototypes_directory = os.path.join(directory, 'prototypes')
        if not os.path.exists(prototypes_directory):
            os.makedirs(prototypes_directory)

        warpings_directory = os.path.join(directory, 'warpings')
        if not os.path.exists(warpings_directory):
            os.makedirs(warpings_directory)

        pois_directory = os.path.join(directory, 'pois')
        if not os.path.exists(pois_directory):
            os.makedirs(pois_directory)

        poi_hists_directory = os.path.join(directory, 'poi_hists')
        if not os.path.exists(poi_hists_directory):
            os.makedirs(poi_hists_directory)

        entropies_file = os.path.join(directory, 'entropies.csv')

        print '> Saving regions to directory {0}'.format(bed_directory)
        for i, c in enumerate(self._clusters):

            if c.regions:
                filename = os.path.join(bed_directory, 'cluster_{0}.bed'.format(i + 1))
                c.save_as_encode_track(filename, track_name='DGWCluster{0}'.format(i + 1),
                                       track_description='DGW Cluster {0}/{1} (cut level: {2})'.format(i + 1,
                                                                                                       len(self._clusters),
                                                                                                       self.cut_level))
            else:
                filename = os.path.join(bed_directory, 'cluster_{0}.txt'.format(i + 1))
                c.save_as_list_of_indices(filename)



        print '> Saving heatmaps to directory {0}'.format(heatmaps_directory)
        for i, c in enumerate(self._clusters):

            filename_reg = os.path.join(heatmaps_directory, 'cluster-{0}.pdf'.format(i + 1))
            f = self._plot_regular_heatmap_on_figure(c)
            f.savefig(filename_reg)
            plt.close(f)


            filename_projected = os.path.join(heatmaps_directory, 'cluster-warped-{0}.pdf'.format(i + 1))
            f = self._plot_projected_heatmap_on_figure(c)
            f.savefig(filename_projected)
            plt.close(f)

        print '> Saving prototypes to directory {0}'.format(prototypes_directory)
        for i, c in enumerate(self._clusters):
            filename_prototype = os.path.join(prototypes_directory, 'cluster-{0}-prototype.pdf'.format(i + 1))
            f = self._plot_prototype_on_figure(c)
            f.savefig(filename_prototype)
            plt.close(f)

            filename_text_prototype = os.path.join(prototypes_directory, 'cluster-{0}-prototype.tsv'.format(i+1))
            c.save_prototype_to_text(filename_text_prototype)

            filename_conservation = os.path.join(prototypes_directory, 'cluster-{0}-conservation.tsv'.format(i+1))
            c.save_conservation_coefficient_as_text(filename_conservation)

        print '> Saving warpings to directory {0}'.format(warpings_directory)
        for i, c in enumerate(self._clusters):
            filename_data = os.path.join(warpings_directory, 'cluster-{0}-warpings.tsv.gz'.format(i + 1))

            c.save_warpings_to_file(filename_data)

            filename_data_conservation = os.path.join(warpings_directory, 'cluster-{0}-warping-conservation.tsv.gz'.format(i+1))
            c.save_warping_conservation_data_to_file(filename_data_conservation)

        print '> Saving POIs to directory {0}'.format(pois_directory)
        for i, c in enumerate(self._clusters):
            filename_data = os.path.join(pois_directory, 'cluster-{0}-pois.tsv.gz'.format(i + 1))
            c.save_pois_to_file(filename_data)

        print '> Saving POI histograms to directory {0}'.format(poi_hists_directory)
        for i, c in enumerate(self._clusters):
            filename_data = os.path.join(poi_hists_directory, 'cluster-{0}-poi-hist'.format(i + 1))
            c.save_poi_histograms_to_file(filename_data)

        print '> Saving entropies to file: {0}'.format(entropies_file)
        self.save_entropies(entropies_file)

        print '> Finished saving clusters'

    def save_entropies(self, file_):
        entropies = []
        for i, c in enumerate(self._clusters):
            cluster_entropies = c.poi_entropies()
            if cluster_entropies is not None:
                cluster_entropies = cluster_entropies
                cluster_entropies['cluster'] = i+1
                cluster_entropies['cluster_size'] = c.n_items

                entropies.append(cluster_entropies)

        if entropies:
            entropies = pd.concat(entropies).reset_index().set_index(['poi_file', 'cluster'])
            entropies.sort_index()
            entropies.to_csv(file_)

    def _callback_save(self, event):
        debug('Save clicked')

        directory = self.output_directory
        self.save_clusters(directory)

    def add_button(self, text, callback):
        buttons_left = 0.03
        buttons_bottom = 1 - 0.1
        button_width = 0.14
        button_height = 0.05
        button_spacing = 0.02

        ax_button = plt.axes([buttons_left + (button_width + button_spacing) * len(self._buttons),
                              buttons_bottom, button_width, button_height])
        self._buttons.append(Button(ax_button, text))
        self._buttons[-1].on_clicked(callback)

    def _plot_prototype_on_figure(self, cluster):
        f = plt.figure()
        plt.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
        cluster.prototype.plot(ax=f.gca())

        return f

    def _plot_regular_heatmap_on_figure(self, cluster):
        f = plt.figure(figsize=(11.7, 8.3))
        plt.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
        self._enlarged_heatmap_axis = plt.gca()

        cluster.data.plot_heatmap(horizontal_grid=True,
                                          sort_by=None, highlighted_points=cluster.points_of_interest, highlight_colours=self.highlight_colours)

        f.canvas.mpl_connect('button_press_event', self._onclick_listener)
        return f

    def _plot_projected_heatmap_on_figure(self, cluster):
        f = plt.figure(figsize=(11.7, 8.3))
        plt.subplots_adjust(left=0.05, bottom=0.05, top=0.95, right=0.95)
        self._enlarged_projected_heatmap_axis = plt.gca()

        cluster.projected_data.plot_heatmap(horizontal_grid=True,
                                                    sort_by=None, highlighted_points=cluster.tracked_points_of_interest,
                                                    highlight_colours=self.highlight_colours)

        f.canvas.mpl_connect('button_press_event', self._onclick_listener)

        return f

    def _enlarge_heatmaps(self, event):
        self._plot_regular_heatmap_on_figure(self.current_cluster())
        plt.show()

        self._plot_projected_heatmap_on_figure(self.current_cluster())
        plt.show()

    def _view_poi_histogram(self, event):
        current_cluster = self.current_cluster()
        if not current_cluster.points_of_interest:
            plt.figure()
            plt.figtext(0.1, 0.1, 'Sorry, no points of interest for this cluster. Try running explorer again with -poi poi_dataset.bed')
            plt.show()
            return

        # Pandas tend to crash if subplots is True when only one dimension
        number_poi_collections = len(current_cluster.points_of_interest_histogram.columns)
        axes = current_cluster.points_of_interest_histogram.plot(kind='bar', subplots=(number_poi_collections > 1))
        if number_poi_collections == 1:
            axes = [axes]
        for ax in axes:
            ax.set_ylabel('Number of regions')
        plt.xlabel('Bin')
        plt.suptitle('Original points of interest')
        axes = current_cluster.tracked_points_histogram.plot(kind='bar', subplots=(number_poi_collections > 1))
        if number_poi_collections == 1:
            axes = [axes]
        for ax in axes:
            ax.set_ylabel('Number of regions')

        plt.xlabel('Warped bin')
        plt.suptitle('Warped points of interest')
        plt.show()

    def _toggle_warping_conservation_view(self, event):
        self._warping_conservation_view = not self._warping_conservation_view
        self.draw()


    def create_buttons(self):
        self.add_button('Previous', self._callback_previous)
        self.add_button('Next', self._callback_next)

        self.add_button('Save clusters', self._callback_save)

        self.add_button('Enlarge heatmaps', self._enlarge_heatmaps)
        self.add_button('View POI Histogram', self._view_poi_histogram)
        self.add_button('Toggle conservation view', self._toggle_warping_conservation_view)

    def _plot_item_in_figure(self, index):
        data = self.current_cluster().data.ix[index]
        projected_data = self.current_cluster().projected_data.ix[index]
        prototype = self.current_cluster().prototype

        poi = self.current_cluster().points_of_interest.get(index, {})
        tracked_poi = self.current_cluster().tracked_points_of_interest.get(index, {})

        plt.figure()
        ax1 = plt.subplot(3, 1, 1)
        data.plot(ax=ax1, legend=False)
        plt.figlegend(*ax1.get_legend_handles_labels(), loc='lower center')
        if poi:
            lim_min, lim_max = ax1.get_ylim()
            height = (lim_max - lim_min) / 20
            points_plotted_on = defaultdict(lambda: 0)
            for key, value in poi.iteritems():
                colour = self.highlight_colours[key]
                for point in value:
                    items_on_current_point = points_plotted_on[point]
                    ax1.add_patch(Rectangle((point, lim_min + (height*items_on_current_point)), width=1, height=height, facecolor=colour, edgecolor='k'))
                    points_plotted_on[point] += 1
        plt.title('Original')

        ax2 = plt.subplot(3, 1, 2, sharey=ax1)
        prototype.plot(ax=ax2, legend=False)
        plt.title('Cluster Prototype')

        ax3 = plt.subplot(3, 1, 3, sharey=ax1)
        projected_data.plot(ax=ax3, legend=False)

        if tracked_poi:
            lim_min, lim_max = ax3.get_ylim()
            height = (lim_max - lim_min) / 20
            points_plotted_on = defaultdict(lambda: 0)
            for key, value in tracked_poi.iteritems():
                colour = self.highlight_colours[key]
                for point in value:
                    items_on_current_point = points_plotted_on[point]
                    ax3.add_patch(Rectangle((point, lim_min + (height*items_on_current_point)), width=1, height=height,
                                            facecolor=colour, edgecolor='k'))
                    points_plotted_on[point] += 1

        plt.title('Warped')
        plt.suptitle(index)

        figure_dtw_mappings = visualise_dtw_mappings(data, prototype, dtw_function=self.dtw_function,
                                                     columns=data.columns, sequence_x_label=index,
                                                     sequence_y_label='Cluster Prototype')



    def _onclick_listener(self, event):
        if not self.ax_heatmap().contains(event)[0] and not self.ax_projected_heatmap().contains(event)[0]:
           if (self._enlarged_heatmap_axis is None or not self._enlarged_heatmap_axis.contains(event)[0]) and \
               (self._enlarged_projected_heatmap_axis is None or not self._enlarged_projected_heatmap_axis.contains(event)[0]):
               return

        if event.ydata is None or event.button != 3:
            return

        ylabel = event.inaxes.yaxis.get_major_formatter()(event.ydata)
        try:
            self._plot_item_in_figure(ylabel)
            plt.show()
        except Exception, e:
            error('Something went wrong when plotting {0!r}, got: {1!r}'.format(ylabel, e))
            raise

    def draw(self):

        ax_prototype = plt.subplot(self.gs_prototype())
        plt.cla()
        plt.title('Prototype')
        current_cluster = self.current_cluster()
        current_cluster.prototype.plot(ax=ax_prototype, legend=False)

        plt.figlegend(*ax_prototype.get_legend_handles_labels(), loc='lower center')

        plt.subplot(self.gs_warping_preservation(), sharex=ax_prototype)
        plt.cla()
        wpc_vector = current_cluster.warping_conservation_vector()
        plt.plot(np.arange(0.5, len(wpc_vector), 1), wpc_vector)


        self.title.set_text('Cluster #{0}/{2} ({1} elements)'.format(self._current_cluster_id + 1,
                                                                     current_cluster.n_items,
                                                                     len(self._clusters)))

        plt.subplot(self.gs_heatmap())
        plt.cla()
        self._ax_heatmap = plt.gca()

        shared_axis = current_cluster.data.plot_heatmap(horizontal_grid=True, subplot_spec=self.gs_heatmap(),
                                                        sort_by=None, highlighted_points=current_cluster.points_of_interest,
                                                        highlight_colours=self.highlight_colours)

        # Projections
        projections = current_cluster.projected_data

        ax_projections = plt.subplot(self.gs_projected_mean())
        plt.cla()
        plt.title('Average DTW projection onto prototype')
        projections.mean().plot(ax=ax_projections, legend=False)

        plt.subplot(self.gs_projected_heatmap())
        plt.cla()
        self._ax_projected_heatmap = plt.gca()

        if self._warping_conservation_view:
            current_cluster.projected_data.plot_heatmap(horizontal_grid=True, subplot_spec=self.gs_projected_heatmap(),
                                                    share_y_axis=shared_axis, sort_by=None,
                                                    highlighted_points=current_cluster.tracked_points_of_interest,
                                                    highlight_colours=self.highlight_colours,
                                                    replace_with_dataframe=current_cluster.warping_conservation_data)
        else:
            current_cluster.projected_data.plot_heatmap(horizontal_grid=True, subplot_spec=self.gs_projected_heatmap(),
                                                                share_y_axis=shared_axis, sort_by=None,
                                                                highlighted_points=current_cluster.tracked_points_of_interest,
                                                                highlight_colours=self.highlight_colours)


        # Finally issue a draw command for the plot
        plt.draw()

    def create_figure(self, figsize=(12, 10), interactive=True):
        self._figure = plt.figure(num=None,
                                  figsize=figsize,
                                  facecolor='w', edgecolor='k')
        self.title = plt.suptitle("") # Create text object for title

        if interactive:
            self.create_buttons()
            self._figure.canvas.mpl_connect('button_press_event', self._onclick_listener)

        self.draw()

    def show(self):
        self.create_figure(interactive=True)
        plt.show()

class HierarchicalClusteringViewer(object):
    _hierarchical_clustering_object = None
    _gs_main = None
    _figure = None
    _ax_dendrogram = None
    _line = None
    _cut_xdata = 0
    _sorted_index = None
    _last_clusters = None
    _last_xdata = None

    def __init__(self, hierarchical_clustering_object, output_directory, configuration_file, highlight_colours=None,
                 cut_xdata=0):
        self._hierarchical_clustering_object = hierarchical_clustering_object
        self.output_directory = output_directory
        self.configuration_file = configuration_file
        self.highlight_colours = highlight_colours
        self._cut_xdata = cut_xdata

    @property
    def clusters(self):
        if self._cut_xdata > 0:
            if self._last_xdata != self._cut_xdata:
                index = self.sorted_index
                assert(index is not None)
                self._last_clusters = self._hierarchical_clustering_object.cut_and_resort(self._cut_xdata, index)
                self._last_xdata = self._cut_xdata
            return self._last_clusters
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

    def create_figure(self, figsize=(12, 10), interactive=True):
        plt.figure(num=None, figsize=figsize, facecolor='w', edgecolor='k')
        self._gs_main = gridspec.GridSpec(2, 2, wspace=0, height_ratios=[1, 15])
        self._figure = plt.gcf()
        self._ax_dendrogram = plt.subplot(self.gs_dendrogram, rasterized=True)

        if interactive:
            self._figure.canvas.mpl_connect('button_press_event', self._onclick_listener)
            self.draw_buttons()

    @property
    def sorted_index(self):
        return self._sorted_index

    def draw_dendrogram(self):
        ax_dendrogram = self.ax_dendrogram
        plt.subplot(ax_dendrogram, rasterized=True)
        hc = self.hierarchical_clustering_object

        dendrogram_dict = hc.dendrogram(orientation='right', get_leaves=True,  distance_sort=True,
                                        color_threshold=self._cut_xdata, ax=ax_dendrogram)
        leaves = dendrogram_dict['leaves']

        plt.setp(plt.gca().get_xticklabels(), rotation='vertical', fontsize=7)
        plt.gca().get_yaxis().set_visible(False)

        if self._sorted_index is None:
            self._sorted_index = hc.data.items[leaves]


    def cluster_previewer(self):
        if not self.clusters:
            return

        pw = ClusterPreviewer(self.clusters, self.output_directory, self.configuration_file, self._cut_xdata,
                              highlight_colours=self.highlight_colours, dtw_function=self.hierarchical_clustering_object.dtw_function)

        return pw

    def _callback_preview(self, event):
        pw = self.cluster_previewer()

        # cluster_previewer() returns None if there are no clusters
        if pw:
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
        self.draw_dendrogram()
        self._draw_cut_line()
        hc = self.hierarchical_clustering_object
        sorted_index = self.sorted_index
        assert(sorted_index is not None)

        DENDROGRAM_SCALE = 10  # scipy.cluster.hierarachy.dendrogram scales all y axis values by tenfold for some reason
        hc.data.plot_heatmap(subplot_spec=self.gs_heatmap, no_y_axis=True, sort_by=sorted_index, share_y_axis=self.ax_dendrogram,
                             scale_y_axis=DENDROGRAM_SCALE, highlighted_points=hc.data.points_of_interest, rasterized=True,
                             highlight_colours=self.highlight_colours)

    def savefig(self, filename, create=True):

        if create:
            self.create_figure(interactive=False)
            self.draw()

        plt.tight_layout()
        plt.savefig(filename)

    def show(self):
        self.create_figure(interactive=True)
        self.draw()
        plt.show()


    def _draw_cut_line(self):
        xdata = self._cut_xdata
        if xdata > 0:
            if self._line is not None:
                self._line.set_xdata(xdata)
            else:
                self._line = self._ax_dendrogram.axvline(linestyle='--', x=xdata, color='k')


    def cut_line(self, xdata):

        self._cut_xdata = xdata
        self.draw_dendrogram()
        self._draw_cut_line()
        self._figure.canvas.draw()


    def _onclick_listener(self, event):
        # Allow only in axis, only the left button and only regions with value
        if event.inaxes != self.ax_dendrogram or event.button not in [2, 3] or event.xdata is None:
            return

        self.cut_line(event.xdata)

