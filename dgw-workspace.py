import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle
import dgw
from dgw.cli import Configuration

import sys
from dgw.cluster import add_path_data
import dgw.cluster.visualisation

configuration_loc = 'dgw_config.pickle' if len(sys.argv) < 2 else sys.argv[1]

def load_from_pickle(filename):
    f = open(filename, 'rb')
    try:
        return pickle.load(f)
    finally:
        f.close()

configuration = load_from_pickle(configuration_loc)
assert(isinstance(configuration, Configuration))

if configuration.parsed_regions_filename:
    regions = load_from_pickle(configuration.parsed_regions_filename)
else:
    regions = None

if configuration.processed_dataset_filename:
    dataset = load_from_pickle(configuration.processed_dataset_filename)
else:
    dataset = load_from_pickle(configuration.dataset_filename)
    if configuration.raw_dataset_filename:
        raw_dataset = load_from_pickle(configuration.raw_dataset_filename)

if configuration.pairwise_distances_filename:
    dm = np.load(configuration.pairwise_distances_filename)
if configuration.linkage_filename:
    linkage = np.load(configuration.linkage_filename)
    prototypes = load_from_pickle(configuration.prototypes_filename)
    warping_paths = load_from_pickle(configuration.warping_paths_filename)


    hc = dgw.cluster.analysis.HierarchicalClustering(dataset, regions, linkage_matrix=linkage, prototypes=prototypes,
                                                     dtw_function=configuration.dtw_function)
    add_path_data(hc.tree_nodes_list, hc.num_obs, warping_paths)
    hcv = dgw.cluster.visualisation.HierarchicalClusteringViewer(hc)