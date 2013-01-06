__author__ = 'saulius'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fastcluster
import scipy.cluster.hierarchy as hierarchy

import genes
import cluster.distance.dtw.std as dtw
import cluster.distance.dtw.fastdtw as fastdtw
import view.heatmap as heatmap

KNOWN_GENES = '../data/knownGenes'
K562_H3K4ME3_REP1 = '../data/interesting/broad/K562/wgEncodeBroadHistoneK562H3k4me3StdAlnRep1.bam'
K562_H3K4ME3_REP2 = '../data/interesting/broad/K562/wgEncodeBroadHistoneK562H3k4me3StdAlnRep2.bam'
