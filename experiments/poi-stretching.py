import argparse
from collections import defaultdict
from dgw.cli import StoreFilenameAction
import cPickle as pickle
import random
import matplotlib.pyplot as plt
from dgw.data.containers import AlignmentsData
from dgw.dtw import uniform_scaling_to_length
from dgw.evaluation.resampling import mutate_sequence
import numpy as np
import pandas as pd


def load_from_pickle(filename):
    f = open(filename)
    try:
        return pickle.load(f)
    finally:
        f.close()
def argument_parser():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('raw_dataset', help='Raw dataset to work on', action=StoreFilenameAction)

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
    dataset.plot_heatmap(highlighted_points=dataset.points_of_interest, sort_by=None)
    plt.show()

def stretch_poi_by_length(values, poi, stretch_by):
    """
    Stretches POI in the values provided to extend the POI region by desired length
    :param values:
    :param poi:
    :param stretch_by:
    :return:
    """
    poi_data = values[poi]
    new_len = len(poi_data) + stretch_by
    new_poi = uniform_scaling_to_length(poi_data, new_len)

    before_poi = values[:poi[0]]
    after_poi = values[poi[-1] + 1:]

    new_data = np.concatenate((before_poi, new_poi, after_poi))
    new_poi = np.asarray(range(poi[0], poi[0] + new_len))
    return new_data, new_poi

def generate_new_dataset(dataset, desired_size, max_stretch_length=20, reverse_randomly=True, mutation_probability=0.01):
    """
    Generates new dataset of desired size by randomly stretching the POI regions of the dataset provided up to
    max_stretch_length

    :param dataset:
    :param desired_size:
    :param max_stretch_length:
    :param reverse_randomly: will reverse some of the regions randomly if set to true
    :param mutation_probability: rougly the percentage of items in the sequence  that will mutate and obtain the value
                                 of some other random index (i.e. 0.01 means that around 1% of items will mutate)

    :return:
    """
    points_of_interest = dataset.points_of_interest
    counts = defaultdict(lambda: 0)

    new_dataset = {}
    new_points_of_interest = {}
    for i in range(desired_size):
        item = random.choice(dataset.items)
        data = dataset.ix[item]
        values = data.values
        poi = points_of_interest[item]

        stretch_length = random.randint(0, max_stretch_length)
        new_data, new_poi = stretch_poi_by_length(values, poi, stretch_length)

        new_i = counts[item] + 1
        counts[item] = new_i

        new_ix = '{0}-{1}'.format(item, new_i)

        if reverse_randomly and random.choice([True, False]):
            new_data = new_data[::-1]
            new_poi = (len(new_data) - 1 - new_poi)[::-1]
            new_ix += 'r'


        if mutation_probability is not None:
            new_data = mutate_sequence(new_data, mutation_probability)

        new_data = pd.DataFrame(new_data, columns=data.columns)
        new_dataset[new_ix] = new_data
        new_points_of_interest[new_ix] = new_poi


    new_dataset = pd.Panel(new_dataset)
    ad = AlignmentsData(new_dataset)
    ad.points_of_interest = new_points_of_interest

    return ad

def save_dataset(dataset, output_file):
    """
    Will save the dataset to output file that is ready to be passed in DGW as --processed-dataset
    :param dataset:
    :param output_file
    :return:
    """
    dataset = dataset.to_log_scale()
    f = open(output_file, 'w')
    try:
        pickle.dump(dataset, f)
    except Exception:
        f.close()

print "> Please run this file interactively: ipython -i exon-strateching.py -- [arguments]"
print "> Look at the source code for some useful functions"

parser = argument_parser()
args = parser.parse_args()

raw_dataset = load_from_pickle(args.raw_dataset)
raw_dataset = raw_dataset.drop_no_pois()
print 'Size of raw dataset: {0}'.format(len(raw_dataset))

# Try something like
# subset = get_random_sample(raw_dataset, 5)
# new_dataset = generate_new_dataset(subset, 100)
# plot_heatmap(new_dataset)
# save_dataset(new_dataset, 'somefile.pickle')
