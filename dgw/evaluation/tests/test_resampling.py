import unittest
import numpy as np
from numpy.testing import assert_array_equal
from dgw.evaluation.resampling import extend_point, shrink_to_a_single_point


class TestExtending(unittest.TestCase):

    def test_extend_point(self):
        a = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        ans = np.array([0, 1, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9, 10])

        assert_array_equal(ans, extend_point(a, 4, 3))

        # multi-dim
        a = np.array([[1, 2], [2, 3], [3, 4]])
        ans = np.array([[1, 2], [2, 3], [2, 3], [3, 4]])

        assert_array_equal(ans, extend_point(a, 1, 2))

    def test_extend_point_left_boundary(self):

        a = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])
        ans = np.array([0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8])

        assert_array_equal(ans, extend_point(a, 0, 4))


    def test_extend_point_right_boundary(self):

        a = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])
        ans = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8])

        assert_array_equal(ans, extend_point(a, 8, 4))

    def test_shrink_to_single_point(self):

        a = np.array([0, 1, 2, 3, 3, 3, 3, 4, 5, 6])
        ans = np.array([0, 1, 2, 3, 4, 5, 6])

        assert_array_equal(ans, shrink_to_a_single_point(a, 3, 4))

        b = np.array([0, 1, 2, 3, 8, 9, 10, 4, 5, 6])
        ans = np.array([0, 1, 2, np.mean([3, 8, 9, 10]), 4, 5, 6])

        assert_array_equal(ans, shrink_to_a_single_point(b, 3, 4))

        # multi-dim
        a = np.array([[1, 2], [2, 3], [2, 3], [3, 4]])
        ans = np.array([[1, 2], [2, 3], [3, 4]])

        assert_array_equal(ans, shrink_to_a_single_point(a, 1, 2))

    def test_shrink_to_single_point_boundary(self):

        a = np.array([0, 1, 2, 3, 4, 5, 6, 6])
        ans = np.array([0, 1, 2, 3, 4, 5, 6])  # Ignore points that go out of bound

        assert_array_equal(ans, shrink_to_a_single_point(a, 6, 4))



