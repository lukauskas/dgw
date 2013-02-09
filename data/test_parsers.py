__author__ = 'saulius'
import unittest
from numpy.testing import *
import numpy as np
import pandas as pd
import parsers

# -- Stub classes to simulate pysam behaviour ---
class StubAlignedRead(object):
    def __init__(self, pos, alen, is_reverse):
        self.pos = pos
        self.alen = alen
        self.is_reverse = is_reverse

    @property
    def aend(self):
        return self.pos + self.alen

class StubSamfile(object):
    __fetch_response = None

    def __init__(self, fetch_response):
        self.__fetch_response = fetch_response

    def fetch(self, chromosome, start, end):
        for item in self.__fetch_response:
            if item.pos < end and item.aend >= start:
                yield item

# -- Tests ------------------------------------------------------------

class TestReadExtendFunctions(unittest.TestCase):
    def test_extend_regular_read(self):
        aligned_read = StubAlignedRead(100, 36, False)
        start, end = parsers.extend_read_to(aligned_read, 100)

        self.assertEqual(100, start)
        self.assertEqual(200, end)

    def test_extend_reverse_read(self):
        aligned_read = StubAlignedRead(164, 36, True)
        start, end = parsers.extend_read_to(aligned_read, 100)

        self.assertEqual(100, start)
        self.assertEqual(200, end)

    def test_invalid_extend_raises_exception(self):
        aligned_read = StubAlignedRead(100, 36, True)

        with self.assertRaises(ValueError):
            parsers.extend_read_to(aligned_read, 10) # 0 < 36

class TestReadCountForRegionReading(unittest.TestCase):

    def setUp(self):

        aligned_reads = [ StubAlignedRead(10, 10, True), # A
                          StubAlignedRead(20, 10, True),  # B
                          StubAlignedRead(5, 10, False),   # C
                          StubAlignedRead(25, 10, False)   # D
                          ]

        #01234567890123456789012345678901234
        #..........AAAAAAAAAABBBBBBBBBB.....
        #.....CCCCCCCCCC..........DDDDDDDDDD

        self.samfile = StubSamfile(aligned_reads)

    def test_non_extended_read_no_binning(self):
        peak_data = parsers.read_samfile_region(self.samfile, 'chr1', 3, 40, resolution=1, extend_to=None)
        correct = np.array([0,0,1,1,1,1,1,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,1,1,1,1,1,0,0,0,0,0])

        assert_array_equal(correct, peak_data)

    def test_non_extended_read_binning(self):

        peak_data = parsers.read_samfile_region(self.samfile, 'chr1', 3, 43, resolution=5, extend_to=None)
        correct = np.array([1,2,2,2,2,2,1,0])

        assert_array_equal(correct, peak_data)

    def test_extended_read_no_binning(self):
        peak_data = parsers.read_samfile_region(self.samfile, 'chr1', 3, 40, resolution=1, extend_to=15)

        #012|3456789012345678901234567890123456789|
        #...|..aaaaaAAAAAAAAAA....................|
        #...|............bbbbbBBBBBBBBBB..........|
        #...|..CCCCCCCCCCccccc.....DDDDDDDDDDddddd|

        correct = np.array([0,0,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,1,1,1,1,1,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1])
        assert_array_equal(correct, peak_data)


    def test_extended_read_binning(self):
        peak_data = parsers.read_samfile_region(self.samfile, 'chr1', 3, 43, resolution=5, extend_to=15)

        #012[34567|89012|34567|89012|34567|89012|34567|89012]
        #...[..aaa|aaAAA|AAAAA|AA...|.....|.....|.....|.....]
        #...[.....|.....|..bbb|bbBBB|BBBBB|BB...|.....|.....]
        #...[..CCC|CCCCC|CCccc|cc...|..DDD|DDDDD|DDddd|dd...]

        correct = np.array([2,2,3,3,2,2,1,1])
        assert_array_equal(correct, peak_data)

    def test_extended_read_no_binning_extended_boundaries(self):

        #0123456789012345[67]8901234567890123456789
        #.....aaaaaAAAAAA[AA]AA....................
        #...............b[bb]bbBBBBBBBBBB..........
        #.....CCCCCCCCCCc[cc]cc.....DDDDDDDDDDddddd

        peak_data = parsers.read_samfile_region(self.samfile, 'chr1', 16, 18, resolution=1, extend_to=15)
        correct = np.array([3,3])
        assert_array_equal(correct, peak_data)

        #012345678901234567890123456789012345678[9012]
        #.....aaaaaAAAAAAAAAA...................[....]
        #...............bbbbbBBBBBBBBBB.........[....]
        #.....CCCCCCCCCCccccc.....DDDDDDDDDDdddd[d...]

        peak_data = parsers.read_samfile_region(self.samfile, 'chr1', 39, 43, resolution=1, extend_to=15)
        correct = np.array([1,0,0,0])
        assert_array_equal(correct, peak_data)

        #01234567890123456789012345678901234567890[123]
        #.....aaaaaAAAAAAAAAA.....................[...]
        #...............bbbbbBBBBBBBBBB...........[...]
        #.....CCCCCCCCCCccccc.....DDDDDDDDDDddddd[....]

        peak_data = parsers.read_samfile_region(self.samfile, 'chr1', 41, 44, resolution=1, extend_to=15)
        correct = np.array([0,0,0])
        assert_array_equal(correct, peak_data)

        #[012345]67890123456789012345678901234567890123
        #[.....a]aaaaAAAAAAAAAA........................
        #[......].........bbbbbBBBBBBBBBB..............
        #[.....C]CCCCCCCCCccccc.....DDDDDDDDDDddddd....

        peak_data = parsers.read_samfile_region(self.samfile, 'chr1', 0, 6, resolution=1, extend_to=15)
        correct = np.array([0,0,0,0,0,2])
        assert_array_equal(correct, peak_data)

        #[01234]567890123456789012345678901234567890123
        #[.....]aaaaaAAAAAAAAAA........................
        #[.....]..........bbbbbBBBBBBBBBB..............
        #[.....]CCCCCCCCCCccccc.....DDDDDDDDDDddddd....

        peak_data = parsers.read_samfile_region(self.samfile, 'chr1', 0, 5, resolution=1, extend_to=15)
        correct = np.array([0,0,0,0,0])
        assert_array_equal(correct, peak_data)

    def test_extended_reads_extended_boundaries_binning(self):

        #012345678901234567[8901234|5678901|2345678]90123
        #.....aaaaaAAAAAAAA[AA.....|.......|.......].....
        #...............bbb[bbBBBBB|BBBBB..|.......].....
        #.....CCCCCCCCCCccc[cc.....|DDDDDDD|DDDdddd]d....

        peak_data = parsers.read_samfile_region(self.samfile, 'chr1',18, 39, resolution=7, extend_to=15)
        correct = np.array([3,2,1])
        assert_array_equal(correct, peak_data)

if __name__ == '__main__':
    unittest.main()