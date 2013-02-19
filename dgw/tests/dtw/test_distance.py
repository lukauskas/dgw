import unittest
from math import sqrt
from scipy.spatial.distance import cosine
import numpy as np
from numpy.testing import *

from dgw.dtw.distance import _strip_nans, dtw_std


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

class TestDTWStd(unittest.TestCase):

    def test_multidimensional_dtw(self):

        a = np.array([[1,2,3, np.nan], [7,8,9,np.nan]]).T
        b = np.array([[10,12,14], [13,15,17]]).T

        # DTW should match the points:
        # (1,7) to (10,13) (distance: sqrt(81 + 36 = 117))
        # (2,8) to (12,15) (distance: sqrt(100 + 49 = 149))
        # (3,9) to (14,17) (distance: sqrt(121 + 64 = 185))

        # Euclidean distance is the one worth testing for, as sqeuclidean will be the same
        # for a.T and b.T as well.
        euclid_distance = sqrt(117) + sqrt(149) + (sqrt(185))
        self.assertAlmostEqual(euclid_distance, dtw_std(a, b, metric='euclidean'))

    def test_multidimensional_dtw_cosine(self):

        a = np.array([[1,2,3, np.nan], [7,8,9,np.nan]]).T
        b = np.array([[10,12,14], [13,15,17]]).T

        # DTW should match the points:
        # (1,7) to (10,13)
        # (2,8) to (12,15)
        # (3,9) to (14,17)

        #
        cosine_distance = cosine(a[0], b[0]) + cosine(a[1], b[1]) + cosine(a[2], b[2])
        self.assertAlmostEqual(cosine_distance, dtw_std(a, b, metric='cosine'))




