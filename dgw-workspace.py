import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import dgw.cluster.analysis

import sys
prefix = 'dgw_' if len(sys.argv) < 2 else sys.argv[1]

# TODO: this should load configuration object up instead

datasets = pd.load('{0}datasets.pd'.format(prefix))
dm = np.load('{0}pairwise_distances.npy'.format(prefix))
regions = pd.load('{0}regions.pd'.format(prefix))

hc = dgw.cluster.analysis.HierarchicalClustering(datasets, dm)

