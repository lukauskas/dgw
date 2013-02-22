import numpy as np
from ...dtw.distance import _strip_nans

class DataFilter(object):
    """
    Abstract class for DataFilter
    """
    def is_valid(self, data_item):
        """
        This is the method that should be over-rided by subclasses

        :param data_item: a raw array of sequences to consider
        :type data_item: np.ndarray
        :return:
        """
        return True

class MinNumberOfReadsFilter(DataFilter):
    _min_number_of_reads = None

    def __init__(self, min_number_of_reads):
        self._min_number_of_reads = min_number_of_reads

    def is_valid(self, data_row):
        no_nans_row = _strip_nans(data_row)
        total_number_of_reads = np.sum(no_nans_row)
        return total_number_of_reads >= self._min_number_of_reads

class HighestPileUpFilter(DataFilter):

    _highest_pileup_threshold = None

    def __init__(self, highest_pileup_threshold):
        self._highest_pileup_threshold = highest_pileup_threshold

    def is_valid(self, data_row):
        no_nans_row = _strip_nans(data_row)
        highest_pileup = np.max(no_nans_row)

        return highest_pileup >= self._highest_pileup_threshold
