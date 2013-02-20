import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
import dgw.cluster.analysis
from dgw.cli import Configuration

import sys
configuration_loc = 'dgw_config.pickle' if len(sys.argv) < 2 else sys.argv[1]

def load_from_pickle(filename):
    f = open(filename, 'rb')
    obj = pickle.load(f)
    f.close()
    return obj

# TODO: this should load configuration object up instead

configuration = load_from_pickle(configuration_loc)
assert(isinstance(configuration, Configuration))

datasets = load_from_pickle(configuration.dataset_filename)
if configuration.raw_dataset_filename:
    raw_datasets = load_from_pickle(configuration.raw_dataset_filename)

try:
    dm = np.load(configuration.pairwise_distances_filename)
except IOError:
    print "Pairwise distances file not found"
    dm = None

regions = load_from_pickle(configuration.parsed_regions_filename)

if dm:
    hc = dgw.cluster.analysis.HierarchicalClustering(datasets, dm, dtw_function=configuration.dtw_function)

