__author__ = 'saulius'
import pandas as pd
import pysam
from logging import debug
import numpy as np

def read_bed(bed_file, resolution=1):
    '''
    Parses the bed file specified.
    If resolution is provided it will automatically adjust the start and end tags as to be able to further
    split the file into bins of resolution width.

    :param bed_file:
    :param resolution:
    :return:
    '''
    peaks = pd.read_csv(bed_file, sep="\t", header=None)

    peaks.columns = ['chromosome', 'start', 'end', 'name', 'score']
    peaks = peaks.set_index('name')

    peaks = clip_to_fit_resolution(peaks, resolution)
    return peaks


def clip_to_fit_resolution(regions, resolution=1):
    '''
    Clips the regions provided to be of correct size so they can be segmented into a number of resolution sized bins
    :param regions:
    :param resolution:
    :return:
    '''
    if resolution == 1:
        return regions

    lens = regions['end'] - regions['start']
    lens.name = 'length'
    regions = regions.join(lens)

    new_regions_data = []

    for ix, row in regions.iterrows():
        remainder = row['length'] % resolution

        if remainder == 0:
            offset_needed = 0
        else:
            offset_needed = resolution - remainder

        add_left = offset_needed / 2
        add_right = offset_needed / 2 + offset_needed % 2

        row['start'] -= add_left

        if row['start'] < 0:
            # Check if we accidentally went sub zero
            add_right += -row['start']
            row['start'] = 0

        row['end']   += add_right

        new_regions_data.append(row[['chromosome', 'start', 'end']])

    new_peaks = pd.DataFrame(new_regions_data, index=regions.index)
    return new_peaks

def read_samfile_region(samfile, chromosome, start, end, resolution=1, extend_to=200):
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
        read_end   = end+extend_to-1
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
        read_end   = end+extend_to-1
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


def read_bam(alignments_file, peaks, resolution=25, extend_to=200):
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
            current_peak_data = read_samfile_region(samfile,
                peak['chromosome'],
                peak['start'],
                peak['end'],
                resolution=resolution,
                extend_to=extend_to)
        except ValueError, e:
            # Most likely caused due to a chromosome not existing in BAM
            debug('Ignoring {0!r} because of {1!r}'.format(peak, e))
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

