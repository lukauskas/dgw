import json
import os
import cPickle as pickle
import dgw
from dgw.cluster import add_path_data
from ..dtw.distance import parametrised_dtw_wrapper
import numpy as np

def load_from_pickle(f):
    if isinstance(f, basestring):
        f = open(f, 'rb')
        close = True
    else:
        close = False
    try:
        return pickle.load(f)
    finally:
        if close:
            f.close()

def strict_load(filename):
    try:
        data = load_from_pickle(filename)
        return data
    except Exception, e:
        raise Exception("Failed to read {0!r}, got {1!r}. "
                        "Make sure that you have all files output from DGW in the same directory as the config file"
        .format(filename, e))

def optional_load(filename):
    try:
        data = load_from_pickle(filename)
        return data
    except Exception:
        return None


class Configuration(object):
    """
    Object that stores the major configuration of DGW console applications
    """
    _ignore_directory = None
    _directory = None
    _prefix = None
    _dtw_kwargs = None
    _prototyping_method = None
    _resolution = None

    DEFAULT_FILENAMES = {
        'pairwise_distances': '{prefix}_pairwise_distances.npy',
        'linkage': '{prefix}_linkage.npy',
        'prototypes': '{prefix}_prototypes.pickle',
        'warping_paths': '{prefix}_warping_paths.pickle',
        'config': '{prefix}_config.dgw',
        'missing_regions': '{prefix}_missing_regions.bed',
        'filtered_regions': '{prefix}_filtered_regions.bed',
        'too_short_regions' : '{prefix}_too_short_regions.bed',
        'too_long_regions' : '{prefix}_too_long_regions.bed',
        'no_poi_regions' : '{prefix}_no_poi_regions.bed',
        'dataset': '{prefix}_dataset.pd',
        'raw_dataset': '{prefix}_dataset_raw.pd',
        'regions': '{prefix}_regions.pd',
    }

    FILENAMES = {}

    def __init__(self, args=None, initial_variables=None):
        if args is None and initial_variables is None:
            raise ValueError('Please provide either args or initial_variables')

        if args:
            self._generate_filenames_from_prefix_and_args(args.prefix, args)
            self._dtw_kwargs = self._dtw_kwargs_from_args(args)
            self._prototyping_method = args.prototyping_method
            self._resolution = args.resolution
            self._use_strand_information = args.use_strand_information
        else:
            self.FILENAMES = initial_variables['FILENAMES']
            self._dtw_kwargs = initial_variables['dtw_kwargs']
            self._prototyping_method = initial_variables['prototyping_method']
            self._resolution = initial_variables['resolution']
            self._use_strand_information = initial_variables['use_strand_information']

    @property
    def directory(self):
        return self._directory

    @directory.setter
    def directory(self, value):
        self._directory = value

    @property
    def resolution(self):
        return self._resolution

    @property
    def use_strand_information(self):
        return self._use_strand_information

    def _generate_filenames_from_prefix_and_args(self, prefix, args):
        """
        Generate default filenames from prefix provided.
        Also parsers the directory

        :param prefix:
        :param args:
        :return:
        """
        self.directory = os.path.dirname(prefix)
        prefix = os.path.basename(prefix)

        self._save_default_filename('config', prefix)

        if not args.blank:
            if args.output_pairwise_distances:
                self._save_default_filename('pairwise_distances', prefix)

            self._save_default_filename('linkage', prefix)
            self._save_default_filename('prototypes', prefix)
            self._save_default_filename('warping_paths', prefix)

        if args.regions:
            self._save_default_filename('regions', prefix)
            self._save_default_filename('missing_regions', prefix)
            self._save_default_filename('filtered_regions', prefix)
            self._save_default_filename('too_short_regions', prefix)

            if args.max_bins:
                self._save_default_filename('too_long_regions', prefix)

        self._save_default_filename('dataset', prefix)
        if args.output_raw_dataset:
            self._save_default_filename('raw_dataset', prefix)

        if args.points_of_interest and args.ignore_no_poi_regions:
            self._save_default_filename('no_poi_regions', prefix)


    def _save_default_filename(self, file_type, prefix):
        self.FILENAMES[file_type] = self.DEFAULT_FILENAMES[file_type].format(prefix=prefix)

    def set_ignore_directory(self, value):
        self._ignore_directory = value

    def _get_filename(self, file_type):
        try:
            return os.path.join(self.directory, self.FILENAMES[file_type])
        except KeyError:
            return None

    @property
    def pairwise_distances_filename(self):
        return self._get_filename('pairwise_distances')

    @property
    def linkage_filename(self):
        return self._get_filename('linkage')

    @property
    def prototypes_filename(self):
        return self._get_filename('prototypes')

    @property
    def warping_paths_filename(self):
        return self._get_filename('warping_paths')

    @property
    def configuration_filename(self):
        return self._get_filename('config')

    @property
    def parsed_regions_filename(self):
        return self._get_filename('regions')

    @property
    def dataset_filename(self):
        return self._get_filename('dataset')

    @property
    def raw_dataset_filename(self):
        return self._get_filename('raw_dataset')

    @property
    def missing_regions_filename(self):
        return self._get_filename('missing_regions')

    @property
    def filtered_regions_filename(self):
        return self._get_filename('filtered_regions')

    @property
    def too_short_regions_filename(self):
        return self._get_filename('too_short_regions')

    @property
    def too_long_regions_filename(self):
        return self._get_filename('too_long_regions')

    @property
    def no_poi_regions_filename(self):
        return self._get_filename('no_poi_regions')

    def _dtw_kwargs_from_args(self, args):
        kw = {'metric': args.metric}
        if args.warping_penalty is not None:
            kw['warping_penalty'] = args.warping_penalty
        if args.slanted_band is not None:
            kw['constraint'] = 'slanted_band'
            kw['k'] = args.slanted_band
        if args.scale:
            kw['scale_first'] = True

        kw['normalise'] = not args.no_length_normalisation
        kw['try_reverse'] = not args.no_reverse
        return kw

    @property
    def dtw_kwargs(self):
        return self._dtw_kwargs

    @property
    def dtw_function(self):
        return parametrised_dtw_wrapper(**self.dtw_kwargs)

    @property
    def prototyping_method(self):
        return self._prototyping_method


    def to_json(self, file):
        data = {'FILENAMES': self.FILENAMES,
                'dtw_kwargs': self.dtw_kwargs,
                'prototyping_method': self._prototyping_method,
                'resolution': self._resolution,
                'use_strand_information' : self._use_strand_information}

        return json.dump(data, file)

    @classmethod
    def from_json(cls, file):
        return Configuration(initial_variables=json.load(file))


    @property
    def blank(self):
        return self.linkage_filename is None

    def load_regions(self):

        return optional_load(self.parsed_regions_filename)

    def load_dataset(self):
        return strict_load(self.dataset_filename)

    def load_linkage(self):
        try:
            linkage = np.load(self.linkage_filename)
        except Exception, e:
            raise IOError('Error reading linkage file {0!r}, got {1!r}', format(linkage, e))
        return linkage

    def load_prototypes(self):
        return strict_load(self.prototypes_filename)

    def load_warping_paths(self):
        return strict_load(self.warping_paths_filename)

    def create_hierarchical_clustering_object(self, dataset=None, regions=None):
        if self.blank:
            raise Exception('Cannot create HierarchicalClustering object from a blank runk')

        if dataset is None:
            dataset = self.load_dataset()

        if regions is None:
            regions = self.load_regions()
        linkage = self.load_linkage()
        prototypes = self.load_prototypes()
        dtw_function = self.dtw_function
        prototyping_method = self.prototyping_method
        warping_paths = self.load_warping_paths()

        hc = dgw.cluster.analysis.HierarchicalClustering(dataset, regions, linkage_matrix=linkage, prototypes=prototypes,
                                                         dtw_function=dtw_function,
                                                         prototyping_method=prototyping_method)
        add_path_data(hc.tree_nodes_list, hc.num_obs, warping_paths)

        return hc





def load_configuration_from_file(configuration_file):
    """
    Loads configuration object from configuration file
    :param configuration_file:
    :return:
    """
    if isinstance(configuration_file, basestring):
        configuration_file = open(configuration_file, 'r')
        should_close = True
    else:
        should_close = False

    try:
        configuration = Configuration.from_json(configuration_file)
    except Exception, e:
        raise Exception('Error opening configuration file provided: {0!r}'.format(e))
    finally:
        if should_close:
            configuration_file.close()

    configuration.directory = os.path.dirname(configuration_file.name)
    return configuration