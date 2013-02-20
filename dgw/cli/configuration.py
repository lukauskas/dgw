from dgw.dtw.distance import dtw_std

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
        return '{0}_pairwise_distances.npy'.format(self.args.prefix)
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
    def missing_regions_filename(self):
        return '{0}_missing_regions.pd'.format(self.args.prefix)

    @property
    def dtw_kwargs(self):
        kw = {'metric': self.args.metric}
        if self.args.slanted_band is not None:
            kw['constraint'] = 'slanted_band'
            kw['k'] = self.args.slanted_band
        return kw

    def dtw_function(self):
        f = lambda x, y: dtw_std(x, y, **self.dtw_kwargs)
        return f

