import unittest
from math import sqrt
from scipy.spatial.distance import cosine
import numpy as np
from numpy.testing import *

from dgw.dtw.distance import dtw_std, warping_conservation_vector
from dgw.dtw.utilities import _strip_nans, reverse_sequence


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

class TestReverse(unittest.TestCase):

    def test_single_dimension(self):
        a = np.array([1,2,3,4,5,6])
        assert_array_equal(np.array([6,5,4,3,2,1]), reverse_sequence(a))

    def test_multi_dimension(self):
        a = np.array([[1, 2], [3,4], [5,6]])
        assert_array_equal(np.array([[5,6], [3,4], [1,2]]), reverse_sequence(a))

    def test_single_dimension_with_nans_appended(self):
        a = np.array([1,2,3,4,np.nan, np.nan])
        assert_array_equal(np.array([4, 3, 2, 1, np.nan, np.nan]), reverse_sequence(a))

    def test_multi_dimension_with_nans_appended(self):
        a = np.array([[1,2],[3,4],[np.nan, np.nan]])
        assert_array_equal(np.array([[3,4],[1,2], [np.nan, np.nan]]), reverse_sequence(a))

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

    def test_slanted_band_dtw(self):

        a = np.array([1, 2, 3, 4, 5, 6, 7])
        b = np.array([1, 1, 1, 1, 4, 4, 4])

        # Slanted band should constrain mapping all elements exactly
        correct_dist = 0 + 1 + 2 ** 2 + 3 ** 2 + 1 + 2 ** 2 + 3 ** 2
        dist = dtw_std(a, b, constraint='slanted_band', k=0)

        self.assertEqual(correct_dist, dist)

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

    def test_reverse_dtw(self):

        a = np.array([1, 2, 3, np.nan])
        b = np.array([3, 2, 1])

        self.assertEqual(0, dtw_std(a, b, try_reverse=True))

        dist, cost, path = dtw_std(a, b, try_reverse=True, dist_only=False)

        assert_array_equal([2, 1, 0], path[0])
        assert_array_equal([0, 1, 2], path[1])

    def test_multidimensional_dtw_reverse(self):

        a = np.array([[1,2,3, np.nan], [7,8,9,np.nan]]).T
        b = np.array([[10,12,14], [13,15,17]]).T

        #Reverse:
        b = reverse_sequence(b)

        # DTW should match the points:
        # (1,7) to (10,13) (distance: sqrt(81 + 36 = 117))
        # (2,8) to (12,15) (distance: sqrt(100 + 49 = 149))
        # (3,9) to (14,17) (distance: sqrt(121 + 64 = 185))

        # Euclidean distance is the one worth testing for, as sqeuclidean will be the same
        # for a.T and b.T as well.
        euclid_distance = sqrt(117) + sqrt(149) + (sqrt(185))
        self.assertAlmostEqual(euclid_distance, dtw_std(a, b, metric='euclidean', try_reverse=True))


class TestWarpingConservationComputation(unittest.TestCase):

    def test_warping_conservation_vector_computed_correctly(self):
        """
        Given the following warping path:
        .......__
        .......|.
        ..____/..
        ..|......
        ./.......
        /........

        the only regions without warping are between
        points  0-2 and 6-7 on the Y axis, therefore the warping conservation vector should be
        0|1|2|3|4|5
         2 2 0 0 0
        """
        path = ([0, 1, 2, 2, 3, 4, 5, 6, 7, 7, 8, 9],
                [0, 1, 2, 3, 3, 3, 3, 3, 4, 5, 5, 5])

        assert_array_equal([2, 2, 0, 0, 0], warping_conservation_vector(path))

    def test_warping_conservation_vector_computed_correctly(self):
        # Same as above, but second sequence is reversed

        path = ([0, 1, 2, 2, 3, 4, 5, 6, 7, 7, 8, 9],
                5 - np.array([0, 1, 2, 3, 3, 3, 3, 3, 4, 5, 5, 5]))

        assert_array_equal([0, 0, 0, 2, 2], warping_conservation_vector(path))





