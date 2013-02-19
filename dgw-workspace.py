import pandas as pd
import numpy as np
import dgw.cluster.analysis

datasets = pd.load('dgw_datasets.pd')
dm = pd.load('dgw_pairwise_distances.npy')

hc = dgw.cluster.analysis.HierarchicalClustering(datasets, dm)

