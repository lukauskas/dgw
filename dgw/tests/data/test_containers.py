import unittest
from numpy.testing import assert_array_equal
import pandas as pd
from dgw.data.containers import Regions
import numpy as np

class TestRegionsClipToResolution(unittest.TestCase):

    def test_clipping_when_only_one_bin_present(self):
        regions = Regions(pd.DataFrame( {'chromosome' : ['chr1', 'chr10'], 'start' : [100, 200], 'end' : [117, 220]} ))

        clipped_df = regions.clip_to_resolution(20)

        self.assertEquals('chr1', clipped_df.ix[0]['chromosome'])
        self.assertEquals(99, clipped_df.ix[0]['start'])
        self.assertEquals(119, clipped_df.ix[0]['end'])

        self.assertEquals('chr10', clipped_df.ix[1]['chromosome'])
        self.assertEquals(200, clipped_df.ix[1]['start'])
        self.assertEquals(220, clipped_df.ix[1]['end'])

    def test_clipping_always_greater_or_equal_than_0(self):
        regions = Regions(pd.DataFrame({'chromosome' : ['chr1'], 'start': [5], 'end': [7]}))

        clipped_df = regions.clip_to_resolution(20)

        self.assertEquals('chr1', clipped_df.ix[0]['chromosome'])
        self.assertEquals(0, clipped_df.ix[0]['start'])
        self.assertEquals(20, clipped_df.ix[0]['end'])

    def test_clipping_multiple_bins(self):
        regions = Regions(pd.DataFrame( {'chromosome' : ['chr1', 'chr10'], 'start' : [100, 200], 'end' : [117, 220]} ))

        clipped_df = regions.clip_to_resolution(5)

        self.assertEquals('chr1', clipped_df.ix[0]['chromosome'])
        self.assertEquals(99, clipped_df.ix[0]['start'])
        self.assertEquals(119, clipped_df.ix[0]['end'])

        self.assertEquals('chr10', clipped_df.ix[1]['chromosome'])
        self.assertEquals(200, clipped_df.ix[1]['start'])
        self.assertEquals(220, clipped_df.ix[1]['end'])

    def test_clipping_res_1(self):

        regions = Regions(pd.DataFrame( {'chromosome' : ['chr1', 'chr10'], 'start' : [100, 200], 'end' : [117, 220]} ))

        clipped_df = regions.clip_to_resolution(1)

        self.assertEquals('chr1', clipped_df.ix[0]['chromosome'])
        self.assertEquals(100, clipped_df.ix[0]['start'])
        self.assertEquals(117, clipped_df.ix[0]['end'])

        self.assertEquals('chr10', clipped_df.ix[1]['chromosome'])
        self.assertEquals(200, clipped_df.ix[1]['start'])
        self.assertEquals(220, clipped_df.ix[1]['end'])

    def test_clipping_keeps_the_same_class(self):

        regions = Regions(pd.DataFrame({'chromosome': ['chr1', 'chr10'], 'start': [100, 200], 'end': [117, 220]}))
        clipped_regions = regions.clip_to_resolution(5)
        self.assertTrue(isinstance(clipped_regions, Regions))

class TestRegionsAsBinsOf(unittest.TestCase):

    def test_raises_value_error_on_mismatch(self):
        parent_regions = Regions(pd.DataFrame({'chromosome': ['chr1'], 'start': [100], 'end':[200]}))

        # Wrong chromosome
        poi_regions = Regions(pd.DataFrame({'chromosome': ['chr2'], 'start': [120], 'end':[180]}))
        self.assertRaises(ValueError, poi_regions.as_bins_of, parent_regions)

        # No overlap 1
        poi_regions2 = Regions(pd.DataFrame({'chromosome': ['chr1'], 'start': [20], 'end':[80]}))
        self.assertRaises(ValueError, poi_regions2.as_bins_of, parent_regions)

        # No overlap 2
        poi_regions3 = Regions(pd.DataFrame({'chromosome': ['chr1'], 'start': [220], 'end':[230]}))
        self.assertRaises(ValueError, poi_regions3.as_bins_of, parent_regions)

    def assert_numpy_dicts_equal(self, correct_dict, other_dict):

        self.assertEqual(correct_dict.keys(), other_dict.keys())

        for key in correct_dict:
            c = correct_dict[key]
            o = other_dict[key]

            assert_array_equal(c, o)

    def test_bins_calculated_correctly(self):
        parent_regions = Regions(pd.DataFrame({'chromosome': ['chr1', 'chr10'], 'start': [100, 200], 'end': [160, 220]}))

        # Edges
        edge_regions = Regions(pd.DataFrame({'chromosome': ['chr1', 'chr10'], 'start': [100, 215], 'end': [110, 220]}))
        correct_ans = {0: np.array([0, 1]), 1: np.array([3])}
        self.assert_numpy_dicts_equal(correct_ans, edge_regions.as_bins_of(parent_regions, resolution=5))

        # Middle
        mid_regions = Regions(pd.DataFrame({'chromosome': ['chr1', 'chr10'], 'start': [140, 208], 'end': [145, 212]}))
        correct_ans = {0: np.array([8]), 1: np.array([1, 2])}
        self.assert_numpy_dicts_equal(correct_ans, mid_regions.as_bins_of(parent_regions, resolution=5))

        # Middle, without one index
        mid_regions2 = Regions(pd.DataFrame({'chromosome': ['chr10'], 'start': [208], 'end': [212]}, index=[1]))
        correct_ans = {1: np.array([1, 2])}
        self.assert_numpy_dicts_equal(correct_ans, mid_regions2.as_bins_of(parent_regions, resolution=5))


        length_one_regions = Regions(pd.DataFrame({'chromosome' : ['chr1'], 'start': [140], 'end': [141]}))
        correct_ans = {0: np.array([8])}
        self.assert_numpy_dicts_equal(correct_ans, length_one_regions.as_bins_of(parent_regions, resolution=5))

    def test_named_indices(self):

        parent_regions = Regions(pd.DataFrame({'chromosome': ['chr1', 'chr10'], 'start': [100, 200], 'end': [160, 220]},
                                 index=['foo', 'bar']))

        # Edges
        edge_regions = Regions(pd.DataFrame({'chromosome': ['chr1', 'chr10'], 'start': [100, 215], 'end': [110, 220]},
                                            index=['foo', 'bar']))
        correct_ans = {'foo': np.array([0, 1]), 'bar': np.array([3])}
        self.assert_numpy_dicts_equal(correct_ans, edge_regions.as_bins_of(parent_regions, resolution=5))


        # Middle
        mid_regions = Regions(pd.DataFrame({'chromosome': ['chr1', 'chr10'], 'start': [140, 208], 'end': [145, 212]},
                                           index=['foo', 'bar']))
        correct_ans = {'foo': np.array([8]), 'bar': np.array([1, 2])}
        self.assert_numpy_dicts_equal(correct_ans, mid_regions.as_bins_of(parent_regions, resolution=5))


        # Middle, without one index
        mid_regions2 = Regions(pd.DataFrame({'chromosome': ['chr10'], 'start': [208], 'end': [212]}, index=['bar']))
        correct_ans = {'bar': np.array([1, 2])}
        self.assert_numpy_dicts_equal(correct_ans, mid_regions2.as_bins_of(parent_regions, resolution=5))


        length_one_regions = Regions(pd.DataFrame({'chromosome': ['chr1'], 'start': [140], 'end': [141]}, index=['foo']))
        correct_ans = {'foo': np.array([8])}
        self.assert_numpy_dicts_equal(correct_ans, length_one_regions.as_bins_of(parent_regions, resolution=5))

        regions_not_in_parent = Regions(pd.DataFrame({'chromosome': ['chr1', 'chr2'], 'start': [140, 200], 'end': [141, 300]},
                                                     index=['foo', 'baz']))
        correct_ans = {'foo': np.array([8])}
        self.assert_numpy_dicts_equal(correct_ans, regions_not_in_parent.as_bins_of(parent_regions, resolution=5))



class TestCommonOperationsDoNotBreakRegions(unittest.TestCase):
    """
    Pretty much tests whether the Issue #2859 in pandas is fixed.
    That is if DataFrame object started playing nicely with subclasses
    """

    def test_slicing_works(self):
        regions = Regions(pd.DataFrame({'chromosome': ['chr1', 'chr10'], 'start': [100, 200], 'end': [117, 220]}))
        self.assertTrue(isinstance(regions.head(1), Regions))
        self.assertTrue(isinstance(regions[:1], Regions))

