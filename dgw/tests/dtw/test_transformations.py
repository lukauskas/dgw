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

        # When scaling from size 7 to size 3
        # each scaled element will consist of 7/3 = 2.333 elements in longer sequence:
        # Thus short[0] is equivalent to average between 0 and 2.333 elements in the long sequence
        #      short[1] is average of elements 2.333 - 4.666 in the long seq
        #      short[2] is average of elements 4.666 - 7 in the long seq

        # Average of fractional elements can be computed by assuming the element "contributes" to average
        # only fractionally e.g. if the interval is [1.2, 3.6) then the avg should be computed as:
        # (1-0.2) * 2 + 1 * 3 + 0.6 * 4 / (0.8 + 1 + 0.6)

        ans = np.array([(1 + 2 + (1 / 3.0) * 3) / (1 + 1 + (1 / 3.0)),  # [0, 2.333)
                        ((1 - (1 / 3.0)) * 3 + 4 + (2 / 3.0) * 5) / (1 - (1 / 3.0) + 1 + (2 / 3.0)),  # [2.333; 4.666]
                        ((1 - (2 / 3.0)) * 5 + 6 + 7) / (1 - (2 / 3.0) + 1 + 1)])  # [4.666; 7)

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

        # When scaling from size 7 to size 3
        # each scaled element will consist of 7/3 = 2.333 elements in longer sequence:
        # Thus short[0] is equivalent to average between 0 and 2.333 elements in the long sequence
        #      short[1] is average of elements 2.333 - 4.666 in the long seq
        #      short[2] is average of elements 4.666 - 7 in the long seq

        # Average of fractional elements can be computed by assuming the element "contributes" to average
        # only fractionally e.g. if the interval is [1.2, 3.6) then the avg should be computed as:
        # (1-0.2) * 2 + 1 * 3 + 0.6 * 4 / (0.8 + 1 + 0.6)

        one_third = 1 / 3.0
        two_thirds = 2 / 3.0
        ans = np.array([[(1 + 2 + one_third * 3) / (1 + 1 + one_third),
                             (8 + 9 + one_third * 10.0) / (1 + 1 + one_third)], # [0, 2.333)
                        [((1 - one_third) * 3 + 4 + two_thirds * 5) / (1 - one_third + 1 + two_thirds),
                             ((1 - one_third) * 10 + 11 + two_thirds * 12) / (1 - one_third + 1 + two_thirds)],  # [2.333; 4.666]
                        [((1 - two_thirds) * 5 + 6 + 7) / (1 - two_thirds + 1 + 1),
                             ((1 - two_thirds) * 12 + 13 + 14) / (1 - two_thirds + 1 + 1)]])  # [4.666; 7)

        result = uniform_shrinking_to_length(b, 3)
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


class TestPathAveraging(unittest.TestCase):

    def test_one_dim(self):
        a = np.array([1, 2, 3])
        b = np.array([5, 8, 9, 11])

        path = (np.array([0, 1, 2, 2]),
                np.array([0, 1, 2, 3]))

        correct_ans = np.array([4, 5, 6, 7])
        average_path = dtw_path_averaging(a, b, path=path)
        assert_array_equal(correct_ans, average_path)

        path2 = (np.array([0, 1, 1, 2]),
                 np.array([0, 1, 2, 3]))

        correct_ans2 = np.array([4, 5, (2 + 9) / 2.0, 7])
        average_path2 = dtw_path_averaging(a, b, path=path2)
        assert_array_equal(correct_ans2, average_path2)

    def test_multi_dim(self):

        a = np.array([[1, 21], [2, 22], [3, 23]])
        b = np.array([[5, 35], [8, 38], [9, 39], [11, 41]])

        path = (np.array([0, 1, 2, 2]),
                np.array([0, 1, 2, 3]))

        correct_ans = np.array([[4, (35 + 21) / 2.0],
                                [5, (38 + 22) / 2.0],
                                [6, ((39 + 23) / 2.0)],
                                [7, (41 + 23) / 2.0]])

        average_path = dtw_path_averaging(a, b, path=path)
        assert_array_equal(correct_ans, average_path)

        path2 = (np.array([0, 1, 1, 2]),
                 np.array([0, 1, 2, 3]))

        correct_ans2 = np.array([[4, (35 + 21) / 2.0],
                                [5, (38 + 22) / 2.0],
                                [(2 + 9) / 2.0, ((39 + 22) / 2.0)],
                                [7, (41 + 23) / 2.0]])
        average_path2 = dtw_path_averaging(a, b, path=path2)
        assert_array_equal(correct_ans2, average_path2)


class TestSdtwAveraging(unittest.TestCase):

    def test_one_dim_weights_set_to_one(self):
        a = np.array([1, 2, 3])
        b = np.array([5, 8, 9, 11])

        path = (np.array([0, 1, 2, 2]),
                np.array([0, 1, 2, 3]))

        # Should be the same as scaled regular path averaging when weights = 1
        correct_ans = uniform_shrinking_to_length(np.array([4, 5, 6, 7]), 4)
        average_path = sdtw_averaging(a, b, 1, 1, path=path)

        assert_array_equal(correct_ans, average_path)

        path2 = (np.array([0, 1, 1, 2]),
                 np.array([0, 1, 2, 3]))

        correct_ans2 = uniform_shrinking_to_length(np.array([4, 5, (2 + 9) / 2.0, 7]), 4)
        average_path2 = sdtw_averaging(a, b, 1, 1, path=path2)
        assert_array_equal(correct_ans2, average_path2)

    def test_multi_dim_weights_set_to_one(self):

        a = np.array([[1, 21], [2, 22], [3, 23]])
        b = np.array([[5, 35], [8, 38], [9, 39], [11, 41]])

        path = (np.array([0, 1, 2, 2]),
                np.array([0, 1, 2, 3]))

        correct_ans = np.array([[4, (35 + 21) / 2.0],
                                [5, (38 + 22) / 2.0],
                                [6, ((39 + 23) / 2.0)],
                                [7, (41 + 23) / 2.0]])

        correct_ans = uniform_shrinking_to_length(correct_ans, 4)
        average_path = sdtw_averaging(a, b, 1, 1, path=path)
        assert_array_equal(correct_ans, average_path)

        path2 = (np.array([0, 1, 1, 2]),
                 np.array([0, 1, 2, 3]))

        correct_ans2 = np.array([[4, (35 + 21) / 2.0],
                                 [5, (38 + 22) / 2.0],
                                 [(2 + 9) / 2.0, ((39 + 22) / 2.0)],
                                 [7, (41 + 23) / 2.0]])
        correct_ans2 = uniform_shrinking_to_length(correct_ans2, 4)
        average_path2 = sdtw_averaging(a, b, 1, 1, path=path2)
        assert_array_equal(correct_ans2, average_path2)

    def test_one_dim_weights_more_than_one(self):
        a = np.array([1, 2, 3])
        b = np.array([5, 8, 9, 11])

        path = (np.array([0, 1, 2, 2]),
                np.array([0, 1, 2, 3]))

        correct_ans = np.array([4, 4, 4, 4, 4,  # Diagonal
                                5, 5, 5, 5, 5,  # Diagonal step
                                6, 6, 6, 6, 6,  # Diagonal step
                                7, 7, 7])       # B moved only

        correct_ans = uniform_shrinking_to_length(correct_ans, 4)
        average_path = sdtw_averaging(a, b, 3, 7, path=path)
        assert_array_equal(correct_ans, average_path)


    def test_multi_dim_weights_set_to_more_than_one(self):

        a = np.array([[1, 21], [2, 22], [3, 23]])
        b = np.array([[5, 35], [8, 38], [9, 39], [11, 41]])

        path = (np.array([0, 1, 2, 2]),
                np.array([0, 1, 2, 3]))

        correct_ans = np.array([[4, (35 + 21) / 2.0], [4, (35 + 21) / 2.0], [4, (35 + 21) / 2.0], [4, (35 + 21) / 2.0],
                                [5, (38 + 22) / 2.0], [5, (38 + 22) / 2.0], [5, (38 + 22) / 2.0], [5, (38 + 22) / 2.0],
                                [6, ((39 + 23) / 2.0)], [6, ((39 + 23) / 2.0)], [6, ((39 + 23) / 2.0)], [6, ((39 + 23) / 2.0)],
                                [7, (41 + 23) / 2.0], [7, (41 + 23) / 2.0], [7, (41 + 23) / 2.0]])

        correct_ans = uniform_shrinking_to_length(correct_ans, 4)
        average_path = sdtw_averaging(a, b, 3, 5, path=path)
        assert_array_equal(correct_ans, average_path)

