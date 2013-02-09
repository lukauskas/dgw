import unittest
import numpy as np
from numpy.testing import *
from distance import _strip_nans

class TestStripNans(unittest.TestCase):

    def test_single_dimension_no_strip_needed(self):
        x = np.array([1,2,3,4,5,6], dtype=float)
        assert_array_equal(x, _strip_nans(x))


    def test_more_dims_no_strip_needed(self):
        x = np.array([[1,2,3,4,5,6], [7,8,9,10,11,12], [13,14,15,16,17,18]], dtype=float).T
        assert_array_equal(x, _strip_nans(x))

    def test_single_dimension_with_strip(self):
        x = np.array([1,2,3,4,5,np.nan, np.nan], dtype=float)
        assert_array_equal(np.array([1,2,3,4,5], dtype=float), _strip_nans(x))

    def test_multi_dimension_with_strip(self):
        x = np.array([[1,2,3,np.nan, np.nan, np.nan], [7,8,9,np.nan, np.nan, np.nan], [13,14,15,np.nan,np.nan,np.nan]], dtype=float).T
        correct = np.array([[1,2,3], [7,8,9], [13, 14, 15]], dtype=float).T
        assert_array_equal(correct, _strip_nans(x))

    def test_multi_dimension_inconsistent_nans(self):
        x = np.array([[1,2,3,4,5,np.nan], [7,8,9,10,11,np.nan], [13,14,15,16,17,18]], dtype=float).T
        self.assertRaises(ValueError, _strip_nans, x)

        x = np.array([[1,2,3,4,5,np.nan], [7,8,9,10,11,np.nan], [13,14,15,16,np.nan,18]], dtype=float).T
        self.assertRaises(ValueError, _strip_nans, x)