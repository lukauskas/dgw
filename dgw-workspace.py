import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import dgw.cluster.analysis

datasets = pd.load('dgw_datasets.pd')
dm = np.load('dgw_pairwise_distances.npy')

hc = dgw.cluster.analysis.HierarchicalClustering(datasets, dm)

