import numpy as np
from pandas.core.nanops import nansum, nanmax

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
        return nansum(data_row) >= self._min_number_of_reads

class HighestPileUpFilter(DataFilter):

    _highest_pileup_threshold = None

    def __init__(self, highest_pileup_threshold):
        self._highest_pileup_threshold = highest_pileup_threshold

    def is_valid(self, data_row):
        return nanmax(data_row) >= self._highest_pileup_threshold
