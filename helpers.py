from math import floor, ceil, factorial
import pandas as pd
import pysam
import numpy as np
import scipy.spatial.distance
from itertools import combinations, izip
from scipy.spatial.distance import num_obs_dm, num_obs_y


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



def read_bed(bed_file, resolution=1):
    peaks = pd.read_csv(bed_file, sep="\t", header=None)

    peaks.columns = ['chromosome', 'start', 'end', 'name', 'score']
    peaks = peaks.set_index('name')

    peaks = clip_to_fit_resolution2(peaks, resolution)
    return peaks

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

def clip_to_fit_resolution(peaks, resolution=1):

    if resolution == 1:
        return peaks

    peaks['start'] = peaks['start'].map(lambda x : resolution * int(floor(float(x) / resolution)) )
    peaks['end']   = peaks['end'].map(lambda x : resolution * int(ceil(float(x) / resolution)) )

    return peaks

def clip_to_fit_resolution2(peaks, resolution=1):
    if resolution == 1:
        return peaks

    lens = peaks['end'] - peaks['start']
    lens.name = 'length'
    peaks = peaks.join(lens)

    new_peaks_data = []
    new_peaks_index = []

    for ix, row in peaks.iterrows():
        remainder = row['length'] % resolution

        if remainder == 0:
            offset_needed = 0
        else:
            offset_needed = resolution - remainder

        add_left = offset_needed / 2
        add_right = offset_needed / 2 + offset_needed % 2


        row['start'] -= add_left
        row['end']   += add_right

        new_peaks_data.append(row[['chromosome', 'start', 'end']])
        new_peaks_index.append(ix)

        assert((row['end'] - row['start']) % resolution == 0)


    new_peaks = pd.DataFrame(new_peaks_data, index=new_peaks_index)
    return new_peaks

def extend_read_to(aligned_read, extend_to):

    if aligned_read.alen > extend_to:
        raise ValueError('AlignedRead {0!r} is already longer than {1}'.format(aligned_read, extend_to))

    if extend_to is None or extend_to == 0:
        alignment_start = aligned_read.pos
        alignment_end   = aligned_read.aend
    else:
        alen = aligned_read.alen
        if alen > extend_to:
            raise ValueError('Alignment length, alen={0} greater than extend_to parameter ({1})'.format(alen, extend_to))
        if not aligned_read.is_reverse:
            alignment_start = aligned_read.pos
            alignment_end   = alignment_start + extend_to
        else:
            alignment_start = aligned_read.aend - extend_to
            alignment_end   = aligned_read.aend

    return alignment_start, alignment_end

def get_read_count_for_region(samfile, chromosome, start, end, resolution=1, extend_to=None):
    '''
        Returns read count data for a samfile.
        Chromosome, start, and end all describe a region of the data.
        This should be zero indexed and the last coordinate is not included, following BED format.

        The read count is returned aggregated by number of base pairs equal to resolution parameter.
        That is if resolution is 10 the data will be returned for every 10 base pairs.

        Note that end-start mod resolution should be zero as the function will not know what to do otherwise.

    :param samfile:
    :param chromosome:
    :param start:
    :param end:
    :param resolution:
    :return:
    '''


    if (end - start) % resolution != 0:
        raise ValueError('The resolution {0} is not valid for peak of length {1}'.format(resolution, end-start))

    # Initialise
    data_len = (end - start) / resolution
    data_buffer = np.zeros(data_len)

    # Fetch all alignments from samfile
    if extend_to is None or extend_to == 0:
        read_start = start
        read_end   = end
    elif extend_to > 0:
        read_start = start-extend_to+1 # +1 because extended reads should overlap with at least one pixel
        read_end   = start+extend_to-1
    else:
        raise ValueError('extend_to should be >= 0')

    alignments = samfile.fetch(chromosome, read_start, read_end)

    for alignment in alignments:
        assert(alignment.aend > alignment.pos)

        if extend_to is not None and extend_to > 0:
            alignment_start, alignment_end = extend_read_to(alignment, extend_to)
            if alignment_end < start or alignment_start >= end:
                continue
        else:
            alignment_start = alignment.pos
            alignment_end   = alignment.aend

        start_bin = max(0, (alignment_start-start) / resolution)
        end_bin   = min(data_len -1, (alignment_end-start - 1) / resolution)

        data_buffer[start_bin:end_bin+1] += 1

    return data_buffer

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

def read_peak_data_from_bam(alignments_file, peaks, resolution=1):
    '''
    Returns data from bam for the specified peaks.
    Peaks should be a pandas.DataFrame object that has 'chromosome', 'start' and 'end' columns.
        Both start and end should be zero-indexed. End coordinate is not included
        This follows BED format of data.

    :param alignments_file:
    :param peaks:
    :param resolution:
    :return:
    '''

    samfile = pysam.Samfile(alignments_file, 'rb')

    peak_data = []
    new_index = []

    assert(isinstance(peaks, pd.DataFrame))
    for index,peak in peaks.iterrows():
        try:
            current_peak_data = get_read_count_for_region(samfile,
                                                        peak['chromosome'],
                                                        peak['start'],
                                                        peak['end'],
                                                        resolution)
        except ValueError:
            # Most likely caused due to a chromosome not existing in BAM
            continue

        peak_data.append(current_peak_data)
        new_index.append(index)

    max_length = max(map(len, peak_data))

    N = len(peak_data)

    # Add NaNs to offsets we do not know peak locations of
    peak_data_np_arr = np.empty([N, max_length])
    for i, x in enumerate(peak_data):
        peak_data_np_arr[i] = np.hstack((x,[np.nan] * (max_length - len(x))))

    peak_data = peak_data_np_arr
#    arr = np.array(peak_data)
    # Create a sparse DataFrame with peak alignments
    sdf = pd.DataFrame(peak_data)
    sdf.index = new_index

    del peak_data

    return sdf



