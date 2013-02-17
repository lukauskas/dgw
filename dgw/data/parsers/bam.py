# coding=utf-8
__author__ = 'saulius'
from logging import debug
import os

import pandas as pd
import pysam
import numpy as np
from dgw.data.containers import AlignmentsData

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


def __read_bam(alignments_filename, regions, resolution=25, extend_to=200):
    """
    Returns data from bam for the regions specified.

    :param alignments_filename: filename of the alignments file to read ss
    :param regions:  `Regions` object defining the regions of interest in the alignments files
    :type regions: Regions
    :param resolution: Resolution at which to read the data
    :type resolution: int
    :param extend_to: The length of alignments will be extended to. Use `None` if no extending needed
    :type extend_to: int or None
    :return: A DataFrame of aggregated reads from BAM file.
    :rtype: pd.DataFrame
    """

    assert(isinstance(regions, Regions))
    regions = regions.clip_to_resolution(resolution)

    source = os.path.basename(alignments_filename)

    samfile = pysam.Samfile(alignments_filename, 'rb')

    peak_data = []
    new_index = []

    assert(isinstance(regions, pd.DataFrame))
    for index,peak in regions.iterrows():
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

def read_bam(alignment_filenames, regions, resolution=25, extend_to=200):
    """
    Reads provided bam files for the data in the specified regions
    :param alignment_filenames: Filenames of the files to read
    :param regions:         Regions to be explored within these files
    :param resolution:      Resolution at which to read the files
    :param extend_to:       Reads will be extended to extend_to base pairs of length, if not None
    :return:
    :rtype: pd.Panel
    """
    # Allow passing either a list of filenames or a single filename
    if type(alignment_filenames) == type('str'):
        alignment_filenames = [alignment_filenames]

    panel_dict = {}

    if len(alignment_filenames) != len(set(alignment_filenames)):
        raise ValueError, 'Some name in alignments_files provided is duplicated. Do you really want to read the file twice?'

    for alignments_file in alignment_filenames:
        name = os.path.basename(alignments_file)
        if name in panel_dict:
            # In case there are files with the same name in differet dirs, use full path
            name = alignments_file

        bam_data = __read_bam(alignments_file, regions, resolution=resolution, extend_to=extend_to)
        panel_dict[name] = bam_data

    panel = pd.Panel(panel_dict)
    # Transpose the panel so the datasets are on the minor axis
    panel.transpose(1, 2, 0)

    return AlignmentsData(panel)
