# coding=utf-8
from itertools import imap
from dgw.data.parsers.filters import MinNumberOfReadsFilter
from dgw.dtw import reverse_sequence

__author__ = 'saulius'
from logging import debug
import os

import pandas as pd
import pysam
import numpy as np
from ..containers import AlignmentsData, Regions

def _alignments_that_fall_within_a_region(samfile, chromosome, start, end, extend_to=200):
    """
    Returns raw alignments that fall within the region of interest in samfile.

    :param samfile:
    :param chromosome:
    :param start:
    :param end:
    :param extend_to:
    :return:
    """
    if start < 0:
        raise ValueError('Read start should be greater than 0, {0} given'.format(start))

    # Fetch all alignments from samfile
    if extend_to is None or extend_to == 0:
        read_start = start
        read_end = end
    elif extend_to > 0:
        read_start = max(start - extend_to + 1, 0) # +1 because extended reads should overlap with at least one pixel
        read_end = end + extend_to-1
    else:
        raise ValueError('extend_to should be >= 0')

    try:
        # Converting to list here as otherwise exception wont be caught
        alignments = list(samfile.fetch(chromosome, read_start, read_end))
    except ValueError:
        if extend_to > 0:
            chromosome_length = samfile.lengths[samfile.references.index(chromosome)]
            if end <= chromosome_length < read_end:  # If the exception was caused by extending
                alignments = samfile.fetch(chromosome, read_start, chromosome_length)
            else:
                raise
        else:
            raise

    return alignments

def _read_samfile_region(samfile, chromosome, start, end, resolution=50, extend_to=200):
    """
        Returns piled up read counts for a samfile.
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
    """
    if (end - start) % resolution != 0:
        raise ValueError('The resolution {0} is not valid for region of length {1}'.format(resolution, end-start))

    # Initialise
    data_len = (end - start) / resolution
    data_buffer = np.zeros(data_len)

    # Fetch all alignments from SamFile
    alignments = _alignments_that_fall_within_a_region(samfile, chromosome, start, end, extend_to=extend_to)
    for alignment in alignments:

        if extend_to is not None and extend_to > 0:
            alignment_start, alignment_end = _extend_read_to(alignment, extend_to)
            if alignment_end < start or alignment_start >= end:
                continue
        else:
            alignment_start = alignment.pos
            alignment_end   = alignment.aend

        start_bin = max(0, (alignment_start-start) / resolution)
        end_bin   = min(data_len -1, (alignment_end-start - 1) / resolution)

        data_buffer[start_bin:end_bin+1] += 1

    return data_buffer

def _extend_read_to(aligned_read, extend_to):

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

def read_bam(alignment_filenames, regions, resolution=50, extend_to=200, data_filters=[MinNumberOfReadsFilter(1)],
            output_removed_indices=False, reverse_negative_strand_regions=False):
    """
    Reads provided bam files for the data in the specified regions
    :param alignment_filenames: Filenames of the files to read
    :param regions:         Regions to be explored within these files
    :param resolution:      Resolution at which to read the files
    :param extend_to:       Reads will be extended to extend_to base pairs of length, if not None
    :param data_filters:    Data filters that will be run on data to determine whether to include the value to list or not
    :param output_removed_indices: function outputs both indices in regions that were not found in dataset, and regions removed by filters
    :param reverse_negative_strand_regions: whether to reverse the regions on negative strand or not
    :return: dataset, [regions_not_in_dataset, regions_removed_by_filter] -- the latter two only if output_removed_indices is set

    :rtype: pd.Panel
    """
    # Allow passing either a list of filenames or a single filename
    if isinstance(alignment_filenames, basestring):
        alignment_filenames = [alignment_filenames]

    if len(alignment_filenames) != len(set(alignment_filenames)):
        raise ValueError('Some name in alignments_files provided is duplicated. '
                         'Do you really want to read the file twice?')

    columns = []
    for alignments_file in alignment_filenames:
        name = os.path.basename(alignments_file)
        if name in columns:
            # In case there are files with the same name in differet dirs, use full path
            name = alignments_file
        columns.append(name)


    samfiles = [pysam.Samfile(alignments_file) for alignments_file in alignment_filenames]

    # Get all references that are in all of the files
    references = reduce(lambda x, y: x & y, [set(s.references) for s in samfiles])

    # Get all regions that are in the dataset by chromosome
    # TODO: make sure the regions removed are somewhere accounted for
    dataset_regions = regions[regions.chromosome.isin(references)]

    indices_not_in_dataset = regions.index - dataset_regions.index

    # Make sure the regions play nicely with the resolution
    dataset_regions = dataset_regions.clip_to_resolution(resolution)

    # find the maximum length of a region (we'll use this to initialise numpy arrays)
    max_len = dataset_regions.lengths.max() / resolution

    dataset = {}

    if reverse_negative_strand_regions:
        if not dataset_regions.has_strand_data:
            raise ValueError('reverse_negative_strand_regions is set to true, yet the regions provided have no strand information.')

    for index, region in dataset_regions.iterrows():
        chromosome = region['chromosome']
        start = region['start']
        end = region['end']
        if reverse_negative_strand_regions:
            strand = region['strand']

        data_arr = np.empty((max_len, len(samfiles)))

        for i, samfile in enumerate(samfiles):
            try:
                region_data = _read_samfile_region(samfile,
                                                   chromosome, start, end,
                                                   resolution=resolution,
                                                   extend_to=extend_to)
            except Exception, e:
                raise IOError('Could not read {0}:{1}-{2} from {4}, got: {3!r}'.format(chromosome, start, end, e, alignment_filenames[i]))

            if reverse_negative_strand_regions and strand == '-':
                region_data = reverse_sequence(region_data)

            padding = [np.nan] * (max_len - len(region_data))
            data_arr[:, i] = np.concatenate((region_data, padding))

        valid = all(imap(lambda f: f.is_valid(data_arr), data_filters))
        if not valid:
            continue

        df = pd.DataFrame(data_arr, columns=columns)
        dataset[index] = df

    panel = pd.Panel(dataset)
    data = AlignmentsData(panel, resolution=resolution)

    filtered_out_indices = dataset_regions.index - data.items

    if output_removed_indices:
        return data, indices_not_in_dataset, filtered_out_indices
    else:
        return data
