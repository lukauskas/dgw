from math import floor, ceil
import pandas as pd
import pysam
import numpy as np

def read_bed(bed_file, resolution=1):
    peaks = pd.read_csv(bed_file, sep="\t", header=None)

    peaks.columns = ['chromosome', 'start', 'end', 'name', 'score']
    peaks = peaks.set_index('name')

    peaks = clip_to_fit_resolution(peaks, resolution)
    return peaks

def join_with_length_information(peaks_df):

    lens = peaks_df['end'] - peaks_df['start']
    lens.name = 'length'

    return peaks_df.join(lens)

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

def get_read_count_for_region(samfile, chromosome, start, end, resolution=1):

    data_arr = []

    # Note this has to be peak.end - 1
    # because samfile.pileup treats ends inclusively
    data_len = (end - start) / resolution
    #    return numpy.random.rand(data_len)

    for i in xrange(data_len):
        start += i * resolution
        end   = start + resolution - 1 # samfile treats ends inclusively
        assert(end > start)

        count = samfile.count(chromosome, start, end)
        data_arr.append(count)

    assert(len(data_arr) > 0)
    assert(len(data_arr) == data_len)
    return data_arr


def read_peak_data_from_bam(alignments_file, peaks, resolution=1):
    '''

    @param alignments_file:
    @param peaks:
    @type peaks: pd.DataFrame
    @param resolution:
    @return:
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
    peak_data = [ x + [np.nan] * (max_length - len(x)) for x in peak_data]

#    arr = np.array(peak_data)
    # Create a sparse DataFrame with peak alignments
    sdf = pd.DataFrame(peak_data)
    sdf.index = new_index

    del peak_data

    return sdf



