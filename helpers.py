from math import floor, ceil, factorial
import pandas as pd
import pysam
import numpy as np
import scipy.spatial.distance
from itertools import combinations, izip

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

def find_prototype(dm, full_index, cluster_indices):

    # If we have only one element in the cluster, it is its own prototype
    if len(cluster_indices) == 1:
        return list(cluster_indices)[0], 0

    sums = {}
    for (i, j), dist in izip(combinations(full_index, 2), dm):
        if i not in cluster_indices:
            continue

        if j not in cluster_indices:
            continue

        try:
            sums[i] += dist
        except KeyError:
            sums[i] = dist

        try:
            sums[j] += dist
        except KeyError:
            sums[j] = dist

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





def get_read_count_for_region(samfile, chromosome, start, end, resolution=1):
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
    alignments = samfile.fetch(chromosome, start, end)

    for alignment in alignments:
        assert(alignment.aend > alignment.pos)
        start_bin = max(0, (alignment.pos-start) / resolution)
        end_bin   = min(data_len -1, (alignment.aend-start - 1) / resolution)

        data_buffer[start_bin:end_bin+1] += 1

    return data_buffer


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



