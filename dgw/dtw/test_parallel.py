import unittest
from scipy.spatial.distance import pdist
import numpy as np
from dgw.dtw.distance import dtw_std
from dgw.dtw.parallel import parallel_pdist
from itertools import combinations
from numpy.testing import assert_array_equal

__author__ = 'saulius'

class TestParallelPdist(unittest.TestCase):

    def setUp(self):
        np.random.seed(42)
        self.sample_data = np.random.randn(10, 20)
        self.sample_data_three_dim = self.sample_data.reshape(10, 20, 1)

    def test_result_returned_same_as_scipy_spatial_distance_pdist_one_process(self):
        correct_ans = pdist(self.sample_data, dtw_std)
        parallel_ans = parallel_pdist(self.sample_data_three_dim, n_processes=1)
        assert_array_equal(correct_ans, parallel_ans)

    def test_result_returned_same_as_scipy_spatial_distance_pdist_one_process_with_metric(self):
        correct_ans = pdist(self.sample_data, lambda x, y: dtw_std(x, y, metric='euclidean'))
        parallel_ans = parallel_pdist(self.sample_data_three_dim, n_processes=1, metric='euclidean')
        assert_array_equal(correct_ans, parallel_ans)

    def test_result_returned_same_as_scipy_spatial_distance_pdist(self):
        correct_ans = pdist(self.sample_data, dtw_std)
        parallel_ans = parallel_pdist(self.sample_data_three_dim)
        assert_array_equal(correct_ans, parallel_ans)

    def test_correct_multidim_result(self):

        data = np.random.randn(10, 16, 2)

        correct_ans = np.array([dtw_std(x, y) for x, y in combinations(data, 2)])

        parallel_ans = parallel_pdist(data, n_processes=1)
        assert_array_equal(correct_ans, parallel_ans)











