import unittest
import StringIO
import cPickle as pickle
import tempfile

import numpy as np

from dgw.data.storage import ChunkedIO


class TestChunkedIOWithFakeFiles(unittest.TestCase):

    def test_chunked_reading(self):
        text = "some sample text"

        f = StringIO.StringIO(text)
        cf = ChunkedIO(f, max_chunk_size=3)  # Specify small chunk size here so it actually uses chunking
        self.assertEqual(text, cf.read())

    def test_chunked_reading_when_size_is_specified(self):

        text = "some sample text"
        f = StringIO.StringIO(text)

        cf = ChunkedIO(f, max_chunk_size=5)

        # Read size smaller than chunk size
        self.assertEqual("so", cf.read(2))
        # Read size equal to chunk size
        self.assertEqual("me sa", cf.read(5))
        # Read size greater than chunk size
        self.assertEqual("mple te", cf.read(7))
        # Read size greater than remaining string
        self.assertEqual("xt", cf.read(10))

    def test_chunked_writing(self):

        text = "some sample text"
        f = StringIO.StringIO()

        cf = ChunkedIO(f, max_chunk_size=5)
        cf.write(text)
        self.assertEqual(text, f.getvalue())

        f.close()
    def test_chunked_writing_smaller_than_chunk_size(self):

        text = "some sample text"
        f = StringIO.StringIO()

        cf = ChunkedIO(f, max_chunk_size=1000)
        cf.write(text)

        self.assertEqual(text, f.getvalue())
        f.close()

    def test_closing(self):

        f = StringIO.StringIO()
        cf = ChunkedIO(f)
        cf.close()

        self.assertTrue(f.closed)

class TestChunkedIOPlaysNiceWithPickle(unittest.TestCase):

    def test_pickle_dump_load_runtrip_works(self):

        o = [1, 2, 3, 4, 5, 6]

        f = StringIO.StringIO()
        cf = ChunkedIO(f, max_chunk_size=2)

        pickle.dump(o, cf, protocol=pickle.HIGHEST_PROTOCOL)
        f.seek(0)

        o2 = pickle.load(cf)

        cf.close()

        self.assertEqual(o, o2)


class TestLargeFilePicklingToDisk(unittest.TestCase):

    def test_large_np_array_pickling(self):
        large_array = np.random.randn(300000000)
        f = tempfile.TemporaryFile()
        cf = ChunkedIO(f)
        pickle.dump(large_array, cf, protocol=pickle.HIGHEST_PROTOCOL)
        f.seek(0)
        large_array2 = pickle.load(cf)
        cf.close()

        self.assertEqual(large_array, large_array2)

    def test_large_np_array_np_save(self):
        large_array = np.random.randn(300000000)
        f = tempfile.TemporaryFile()
        cf = ChunkedIO(f)
        np.save(cf, large_array)
        f.seek(-1)
        large_array2 = pickle.load(cf)
        cf.close()

        self.assertEqual(large_array, large_array2)

