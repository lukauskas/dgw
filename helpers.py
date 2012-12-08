from math import floor, ceil
import pandas as pd
import pysam
import numpy as np

def read_bed(bed_file):
    peaks = pd.read_csv(bed_file, sep="\t", header=None)

    peaks.columns = ['chromosome', 'start', 'end', 'name', 'score']
    peaks = peaks.set_index('name')

    return peaks

def join_with_length_information(peaks_df):

    lens = peaks_df['end'] - peaks_df['start']
    lens.name = 'length'

    return peaks_df.join(lens)


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


def get_peak_data(alignments_file, peaks, resolution=1):
    '''

    @param alignments_file:
    @param peaks:
    @type peaks: pd.DataFrame
    @param resolution:
    @return:
    '''

    samfile = pysam.Samfile(alignments_file, 'rb')

    peak_data = []
    assert(isinstance(peaks, pd.DataFrame))
    for _,peak in peaks.iterrows():
        current_peak_data = get_read_count_for_region(samfile,
                                                    peak['chromosome'],
                                                    peak['start'],
                                                    peak['end'],
                                                    resolution)
        peak_data.append(current_peak_data)

    max_length = max(map(len, peak_data))

    N = len(peak_data)

    # Add NaNs to offsets we do not know peak locations of
    peak_data = [ x + [np.nan] * (max_length - len(x)) for x in peak_data]

    arr = np.array(peak_data)
    # Create a sparse DataFrame with peak alignments
    #sdf = pd.DataFrame(peak_data)
    #sdf.index = peaks.index

    del peak_data

    return arr



