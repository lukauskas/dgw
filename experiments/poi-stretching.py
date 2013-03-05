import dgw
import argparse
from dgw.cli import StoreFilenameAction
import cPickle as pickle
import random
import matplotlib.pyplot as plt


def load_from_pickle(filename):
    f = open(filename)
    try:
        return pickle.load(f)
    finally:
        f.close()
def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('raw_dataset', help='Raw dataset to work on', action=StoreFilenameAction)
    parser.add_argument('-O', '--output', help='Output file', type=argparse.FileType('w'), required=True)

    return parser

def get_random_sample(dataset, n):
    """
    Returns a random sample of the dataset of size n

    :param dataset:
    :param n:
    :return:
    """
    return dataset.ix[random.sample(dataset.items, n)]

def plot_heatmap(dataset):
    """
    Plots heatmap of dataset with POIs superimposed on it
    :param dataset:
    :return:
    """
    dataset.plot_heatmap(highlighted_points=dataset.points_of_interest)
    plt.show()

print "> Please run this file interactively: ipython -i exon-strateching.py -- [arguments]"
print "> Look at the source code for some useful functions"

parser = argument_parser()
args = parser.parse_args()

raw_dataset = load_from_pickle(args.raw_dataset)
raw_dataset = raw_dataset.drop_no_pois()
print 'Size of raw dataset: {0}'.format(len(raw_dataset))

