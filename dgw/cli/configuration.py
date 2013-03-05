import json
import os
from ..dtw.distance import parametrised_dtw_wrapper

class Configuration(object):
    """
    Object that stores the major configuration of DGW console applications
    """
    _ignore_directory = None
    _directory = None
    _prefix = None
    _dtw_kwargs = None
    _prototyping_method = None

    DEFAULT_FILENAMES = {
        'pairwise_distances': '{prefix}_pairwise_distances.npy',
        'linkage': '{prefix}_linkage.npy',
        'prototypes': '{prefix}_prototypes.pickle',
        'warping_paths': '{prefix}_warping_paths.pickle',
        'config': '{prefix}_config.dgw',
        'missing_regions': '{prefix}_missing_regions.pickle',
        'filtered_regions': '{prefix}_filtered_regions.pickle',
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
        else:
            self.FILENAMES = initial_variables['FILENAMES']
            self._dtw_kwargs = initial_variables['dtw_kwargs']
            self._prototyping_method = initial_variables['prototyping_method']

    @property
    def directory(self):
        return self._directory

    @directory.setter
    def directory(self, value):
        self._directory = value


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
        if args.slanted_band is not None:
            kw['constraint'] = 'slanted_band'
            kw['k'] = args.slanted_band
        if args.scale:
            kw['scale_first'] = True

        kw['normalise'] = args.normalise
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
                'prototyping_method': self._prototyping_method}

        return json.dump(data, file)

    @classmethod
    def from_json(cls, file):
        return Configuration(initial_variables=json.load(file))