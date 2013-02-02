from math import floor, ceil, factorial
import pandas as pd
import pysam
import numpy as np
import scipy.spatial.distance
from itertools import combinations, izip
from scipy.spatial.distance import num_obs_dm, num_obs_y

from logging import debug, warn
from data.parsers import read_samfile_region as get_read_count_for_region

def get_read_counts_distribution(regions, samfile):
    def safe_samfile_count(*args, **kwargs):
        try:
            return samfile.count(*args, **kwargs)
        except ValueError:
            return np.nan

    read_counts = map(lambda x : safe_samfile_count(x[1]['chromosome'], x[1]['start'], x[1]['end']), regions.iterrows())

    return pd.Series(read_counts, index=regions.index, name='read_count')

def get_mean_mark(regions, samfile, resolution=1):

    cumulative_mark = None
    count = 0

    for index, data in regions.iterrows():
        read_count = get_read_count_for_region(samfile, data['chromosome'], data['start'], data['end'])

        try:
            cumulative_mark += read_count
        except TypeError:
            # Will happend for the first mark
            cumulative_mark = np.zeros(len(read_count), dtype=float) # Explicity doing np zeros here to make it float
            cumulative_mark += read_count

        count += 1

    if cumulative_mark is not None:
        return cumulative_mark / count
    else:
        return None

def get_highest_stack_distribution(regions, samfile, resolution=1):

    def get_highest_stack_for_item(iter_item):
        index,data = iter_item

        return get_read_count_for_region(samfile, data['chromosome'], data['start'], data['end']).max()

    return pd.Series(map(get_highest_stack_for_item, regions.iterrows()), index=regions.index, name='highest_stack')

def get_variance_of_stack_heights_distribution(regions, samfile, resolution=1):

    def get_variance_of_stack_heights_for_item(iter_item):
        index,data = iter_item

        return get_read_count_for_region(samfile, data['chromosome'], data['start'], data['end']).var()

    return pd.Series(map(get_variance_of_stack_heights_for_item, regions.iterrows()), index=regions.index, name='highest_stack')



def join_with_length_information(peaks_df):

    lens = peaks_df['end'] - peaks_df['start']
    lens.name = 'length'

    return peaks_df.join(lens)

def choose(n, k):
    return factorial(n) / ((factorial(n-k)) * factorial(k))

def read_dm_for_index(dm, index, n=None):

    if n is None:
       n = num_obs_y(dm)

    if index >= n:
        raise ValueError('No index {0} in dm of size {1}'.format(index, n))

    obs = []
    start = 0
    for i in range(index + 1):

        # Each i will have to be compared with n_compared to items:
        n_compared_to = n - i - 1

        if i < index:
            # the comparisons will be in this order
            # (i, i+1), (i, i+2), (i, i+3)
            # So if index == i+1 then offset = 0
            # If index == i+2, offset = 1
            # So on, so offset = index - i - 1
            index_offset = index - i-1
            obs.append(dm[start + index_offset])

            start += n_compared_to
        elif i == index:
            obs.extend(dm[start:start+n_compared_to])

    return np.array(obs)

def find_prototype(dm, full_index, cluster_indices):

    # If we have only one element in the cluster, it is its own prototype
    if len(cluster_indices) == 1:
        return list(cluster_indices)[0], 0

    sums = {}
    full_index_list = list(full_index)
    cluster_indices_is = [ full_index_list.index(x) for x in cluster_indices ]

    for i, ix in zip(cluster_indices_is, cluster_indices):
        sums[ix] = read_dm_for_index(dm, i).sum()

    min_key = None
    min_value = None
    for key, value in sums.iteritems():
        if min_value is None or min_value > value:
            min_value = value
            min_key = key

    return min_key, min_value

def distance_between(dm, n, i, j):
    if i == j:
        return 0
    return dm[len(dm) - choose(n-i,2) + j - i - 1]

def distances_to_previous_peaks(peaks_df):

    distances = {}
    prev_row = None
    for ix, data in peaks_df.iterrows():
        if prev_row is not None and prev_row['chromosome'] == data['chromosome']:
            distances[ix] = data['start'] - prev_row['end']
        else:
            distances[ix] = np.nan

        prev_row = data

    df = pd.Series(distances)
    df.name = 'distance_to_prev_peak'

    del distances
    return df

def join_peaks_that_are_close(peaks_df, proximity_threshold):
    '''
        Joins peaks that are closer than proximity_threshold together.

    @param peaks_df:
    @param proximity_threshold:
    @return:
    '''

    new_peaks_df_data = []
    new_peaks_df_index = []

    new_peak = None
    new_peaks_idx = []

    peaks_df = peaks_df.join(distances_to_previous_peaks(peaks_df))
    for ix, current_peak in peaks_df.iterrows():
        if new_peak is not None \
           and new_peak['chromosome'] == current_peak['chromosome'] \
           and current_peak['distance_to_prev_peak'] <= proximity_threshold:

           # If we can join the peaks, make sure to do so
           new_peak['end'] = current_peak['end']
           new_peaks_idx.append(ix)

        else:
            if new_peak is not None:
                new_peaks_df_data.append(new_peak)
                new_peaks_df_index.append('-'.join(new_peaks_idx))


            new_peak = current_peak[['chromosome', 'start', 'end']]
            new_peaks_idx = [ix]

    if new_peak is not None:
        new_peaks_df_data.append(new_peak)
        new_peaks_df_index.append('-'.join(new_peaks_idx))
        new_peak = None
        new_peaks_idx = []

    return pd.DataFrame(new_peaks_df_data, index=new_peaks_df_index)

def compare_distance_matrices(dm1, dm2):

    def get_closer_distances(df_dm):
        df_dm = df_dm.sort(columns=0)

        closer_distances = {}
        prev_items = set()
        curr_items = set()
        prev_value = None
        for index, value in df_dm.itertuples():
            if prev_value == value:
                closer_distances[index] = prev_items
            elif prev_value < value or prev_value is None:
               prev_items = curr_items.copy()
               prev_value = value
               closer_distances[index] = prev_items
            else:
                raise Exception('Prev_value > value?? prev:{0}, curr:{1}'.format(prev_value, value))

            curr_items.add(index)

        return closer_distances

    # Convert the distance matrices into pd.DataFrame for easier keeping track of indices
    df_dm1 = pd.DataFrame(dm1)
    df_dm2 = pd.DataFrame(dm2)

    closer_distances1 = get_closer_distances(df_dm1)
    closer_distances2 = get_closer_distances(df_dm2)

    overlaps = []
    for i in df_dm1.index:
        overlaps.append(len(closer_distances1[i] & closer_distances2[i]))

    return pd.DataFrame(overlaps, index=df_dm1.index)


