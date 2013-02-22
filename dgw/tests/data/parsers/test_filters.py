import unittest
import numpy as np
from dgw.data.parsers.filters import MinNumberOfReadsFilter, HighestPileUpFilter


class TestHighestPileUpFilter(unittest.TestCase):

    def test_all_zeros(self):
        filter = HighestPileUpFilter(1)
        a = np.array([[0, 0], [0, 0], [0, 0], [np.nan, np.nan]])
        self.assertFalse(filter.is_valid(a))

    def test_equality_passes(self):

        filter = HighestPileUpFilter(1)
        a = np.array([1, 1, 1, 1])
        self.assertTrue(filter.is_valid(a))

    def test_max_value_is_considered_only(self):
        filter = HighestPileUpFilter(5)
        a = np.array([1, 2, 3, 4, 17, 0, 0])
        self.assertTrue(filter.is_valid(a))

        b = np.array([1, 2, 3, 4, 0, 0])
        self.assertFalse(filter.is_valid(b))

    def test_any_dimension_is_fine(self):

        filter = HighestPileUpFilter(5)
        a = np.array([[1, 0], [2, 1], [0, 7], [0, 0], [np.nan, np.nan]])
        self.assertTrue(filter.is_valid(a))

class TestMinNumberOfReadsFilter(unittest.TestCase):

    def test_general_cases(self):
        a = np.array([1, 2, 3, np.nan])
        b = np.array([1, 0, 2, np.nan])
        c = np.array([1, 0, 0, np.nan])

        filter = MinNumberOfReadsFilter(3)
        self.assertTrue(filter.is_valid(a))
        self.assertTrue(filter.is_valid(b))  # Equality should pass
        self.assertFalse(filter.is_valid(c))

    def test_dimensions_are_summed(self):
        a = np.array([[1, 2], [2, 1]])
        filter = MinNumberOfReadsFilter(4)

        self.assertTrue(filter.is_valid(a))

