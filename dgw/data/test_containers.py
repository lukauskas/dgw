import unittest
import pandas as pd
from dgw.data.containers import Regions

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


class TestCommonOperationsDoNotBreakRegions(unittest.TestCase):
    """
    Pretty much tests whether the Issue #2859 in pandas is fixed.
    That is if DataFrame object started playing nicely with subclasses
    """

    def test_slicing_works(self):
        regions = Regions(pd.DataFrame({'chromosome': ['chr1', 'chr10'], 'start': [100, 200], 'end': [117, 220]}))
        self.assertTrue(isinstance(regions.head(1), Regions))
        self.assertTrue(isinstance(regions[:1], Regions))

