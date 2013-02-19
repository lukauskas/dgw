import unittest
import numpy as np
from dgw.dtw.transformations import *
from numpy.testing import assert_array_equal, assert_array_almost_equal


class TestScaling(unittest.TestCase):

    def test_uniform_extending_to_length_one_dim(self):
        a = np.array([1, 2, 3])

        # Easy cases - new dimension is a multiple of previous
        assert_array_equal(np.array([1, 1, 2, 2, 3, 3]), uniform_scaling_to_length(a, 6))

        assert_array_equal(np.array([1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]), uniform_scaling_to_length(a, 12))

        # Harder cases - new dimension is not a multiple of previous
        assert_array_equal(np.array([1, 1, 2, 3]), uniform_scaling_to_length(a, 4))
        assert_array_equal(np.array([1, 1, 2, 2, 3]), uniform_scaling_to_length(a, 5))

    def test_uniform_extending_to_length_two_dim(self):

        a = np.array([[1, 4], [2, 5], [3, 6]])

        # Easy cases -- new dimension is a multiple of previous
        assert_array_equal(np.array([[1, 4], [1, 4], [2, 5], [2, 5], [3, 6], [3, 6]]),
                           uniform_scaling_to_length(a, 6))

        assert_array_equal(np.array([[1, 4], [1, 4], [1, 4], [1, 4],
                                     [2, 5], [2, 5], [2, 5], [2, 5],
                                     [3, 6], [3, 6], [3, 6], [3, 6]]),
                           uniform_scaling_to_length(a, 12))

        # Harder ases - new dimension is not a multiple of previous
        assert_array_equal(np.array([[1, 4], [1, 4], [2, 5], [3, 6]]), uniform_scaling_to_length(a, 4))
        assert_array_equal(np.array([[1, 4], [1, 4], [2, 5], [2, 5], [3, 6]]), uniform_scaling_to_length(a, 5))
        
    def test_uniform_extending_with_nans(self):
        a = np.array([1, 2, 3, np.nan])
        assert_array_equal(np.array([1, 1, 2, 2, 3, 3]), uniform_scaling_to_length(a, 6))

        b = np.array([[1, 4], [2, 5], [3, 6], [np.nan, np.nan]])
        assert_array_equal(np.array([[1, 4], [1, 4], [2, 5], [2, 5], [3, 6], [3, 6]]),
                           uniform_scaling_to_length(b, 6))

    def test_uniform_extending_edge_case_one_element(self):
        a = np.array([1])
        assert_array_equal(np.array([1, 1, 1, 1]), uniform_scaling_to_length(a, 4))

        a2 = np.array([1, np.nan])
        assert_array_equal(np.array([1, 1, 1, 1]), uniform_scaling_to_length(a2, 4))

        b = np.array([[1, 17]])
        assert_array_equal(np.array([[1, 17], [1, 17], [1, 17], [1, 17]]), uniform_scaling_to_length(b, 4))

        b2 = np.array([[1, 17], [np.nan, np.nan]])
        assert_array_equal(np.array([[1, 17], [1, 17], [1, 17], [1, 17]]), uniform_scaling_to_length(b2, 4))

    def test_uniform_extending_edge_case_no_elements(self):
        a = np.array([])
        self.assertRaises(ValueError, uniform_scaling_to_length, a, 4)

        b = np.array([np.nan])
        self.assertRaises(ValueError, uniform_scaling_to_length, b, 4)

    def test_uniform_extending_fails_when_shrinking_asked_for(self):
        a = np.array([1, 2, 3])
        self.assertRaises(ValueError, uniform_scaling_to_length, a, 2)

    def test_uniform_extending_does_nothing_when_length_is_the_same(self):
        a = np.array([1, 2, 3])
        assert_array_equal(a, uniform_scaling_to_length(a, len(a)))

        b = np.array([[1, 4], [2, 5], [3, 6]])
        assert_array_equal(b, uniform_scaling_to_length(b, len(b)))

    def test_uniform_shrinking_to_length_one_dim(self):
        a = np.array([1, 2, 3, 4, 5, 6])

        # Easy cases each item is an average of the blocks
        assert_array_equal(np.array([1.5, 3.5, 5.5]), uniform_shrinking_to_length(a, 3))
        assert_array_equal(np.array([2, 5]), uniform_shrinking_to_length(a, 2))
        assert_array_equal(np.array([3.5]), uniform_shrinking_to_length(a, 1))

        # Harder cases
        b = np.array([1, 2, 3, 4, 5, 6, 7])

        # When scaling from size 6 to size 5
        # each scaled element will consist of 7/3 = 2.333 elements in longer sequence:
        # Thus short[0] is equivalent to average between 0 and 2.333 elements in the long sequence
        #      short[1] is average of elements 2.333 - 4.666 in the long seq
        #      short[2] is average of elements 4.666 - 7 in the long seq

        # Average of fractional elements can be computed by assuming the element "contributes" to average
        # only fractionally e.g. if the interval is [1.2, 3.6) then the avg should be computed as:
        # (1-0.2) * 2 + 1 * 3 + 0.6 * 4 / (0.8 + 1 + 0.6)

        ans = np.array([(1 + 2 + 0.333 * 3) / (1 + 1 + 0.333),  # [0, 2.333)
                        ((1 - 0.333) * 3 + 4 + 0.666 * 5) / (1 - 0.333 + 1 + 0.666),  # [2.333; 4.666]
                        ((1 - 0.666) * 5 + 6 + 7) / (1 - 0.666 + 1 + 1)])  # [4.666; 7)

        result = uniform_shrinking_to_length(b, 3)

        assert_array_almost_equal(ans, result)


    def test_uniform_shrinking_to_length_two_dim(self):
        a = np.array([[1, 7], [2, 8], [3, 9], [4, 10], [5, 11], [6, 12]])

        # Easy cases each item is an average of the blocks
        assert_array_equal(np.array([[1.5, 7.5], [3.5, 9.5], [5.5, 11.5]]), uniform_shrinking_to_length(a, 3))
        assert_array_equal(np.array([[2, 8], [5, 11]]), uniform_shrinking_to_length(a, 2))
        assert_array_equal(np.array([[3.5, 9.5]]), uniform_shrinking_to_length(a, 1))

        # Harder cases
        b = np.array([[1, 8], [2, 9], [3, 10], [4, 11], [5, 12], [6, 13], [7, 14]])

        # When scaling from size 6 to size 5
        # each scaled element will consist of 7/3 = 2.333 elements in longer sequence:
        # Thus short[0] is equivalent to average between 0 and 2.333 elements in the long sequence
        #      short[1] is average of elements 2.333 - 4.666 in the long seq
        #      short[2] is average of elements 4.666 - 7 in the long seq

        # Average of fractional elements can be computed by assuming the element "contributes" to average
        # only fractionally e.g. if the interval is [1.2, 3.6) then the avg should be computed as:
        # (1-0.2) * 2 + 1 * 3 + 0.6 * 4 / (0.8 + 1 + 0.6)

        ans = np.array([[(1 + 2 + 0.333 * 3) / (1 + 1 + 0.333), (8 + 9 + 0.333 * 10) / (1 + 1 + 0.333)]  # [0, 2.333)
                        [((1 - 0.333) * 3 + 4 + 0.666 * 5) / (1 - 0.333 + 1 + 0.666), ((1 - 0.333) * 10 + 11 + 0.666 * 12) / (1 - 0.333 + 1 + 0.666)],  # [2.333; 4.666]
                        [((1 - 0.666) * 5 + 6 + 7) / (1 - 0.666 + 1 + 1), ((1 - 0.666) * 12 + 13 + 14) / (1 - 0.666 + 1 + 1)]])  # [4.666; 7)

        result = uniform_scaling_to_length(b, 3)
        assert_array_almost_equal(ans, result)

    def test_uniform_shrinking_to_length_raises_exception_when_extension_needed(self):

        self.assertRaises(ValueError, uniform_shrinking_to_length, np.array([1, 2, 3]), 6)
        self.assertRaises(ValueError, uniform_shrinking_to_length, np.array([1, 2, 3, np.nan, np.nan]), 4)
        self.assertRaises(ValueError, uniform_shrinking_to_length, np.array([np.nan, np.nan, np.nan, np.nan]), 3)

    def test_uniform_shrinking_does_nothing_when_lengths_equal(self):
        a = np.array([1, 2, 3])
        assert_array_equal(a, uniform_shrinking_to_length(a, 3))

        b = np.array([1, 2, 3, np.nan, np.nan])
        assert_array_equal(a, uniform_shrinking_to_length(b, 3))

    def test_uniform_shrinking_raises_exception_on_invalid_length(self):
        self.assertRaises(ValueError, uniform_shrinking_to_length, np.array([1, 2, 3]), 0)

