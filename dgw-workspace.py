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
    try:
        return pickle.load(f)
    finally:
        f.close()

configuration = load_from_pickle(configuration_loc)
assert(isinstance(configuration, Configuration))

regions = load_from_pickle(configuration.parsed_regions_filename)
dataset = load_from_pickle(configuration.dataset_filename)
if configuration.raw_dataset_filename:
    raw_dataset = load_from_pickle(configuration.raw_dataset_filename)

if configuration.pairwise_distances_filename:
    dm = np.load(configuration.pairwise_distances_filename)
    hc = dgw.cluster.analysis.HierarchicalClustering(dataset, dm, dtw_function=configuration.dtw_function)

