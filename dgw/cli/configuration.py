from dgw.dtw.distance import dtw_std, parametrised_dtw_wrapper


class Configuration(object):
    """
    Object that stores the major configuration of DGW console applications
    """
    _args = None

    def __init__(self, args):
        self._args = args

    @property
    def args(self):
        return self._args

    @property
    def pairwise_distances_filename(self):
        if not self.args.blank:
            return '{0}_pairwise_distances.npy'.format(self.args.prefix)
        else:
            return None

    @property
    def configuration_filename(self):
        return '{0}_config.pickle'.format(self.args.prefix)

    @property
    def parsed_regions_filename(self):
        return '{0}_regions.pd'.format(self.args.prefix)

    @property
    def dataset_filename(self):
        return '{0}_datasets.pd'.format(self.args.prefix)

    @property
    def raw_dataset_filename(self):
        if self.args.output_raw_dataset:
            return '{0}_datasets_raw.pd'.format(self.args.prefix)
        else:
            return None

    @property
    def missing_regions_filename(self):
        return '{0}_missing_regions.pd'.format(self.args.prefix)

    @property
    def dtw_kwargs(self):
        kw = {'metric': self.args.metric}
        if self.args.slanted_band is not None:
            kw['constraint'] = 'slanted_band'
            kw['k'] = self.args.slanted_band
        return kw

    @property
    def dtw_function(self):
        return parametrised_dtw_wrapper((), self.dtw_kwargs)
