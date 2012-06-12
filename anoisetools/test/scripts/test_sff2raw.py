import os
import os.path
import unittest
import tempfile

from anoisetools.scripts import sff2raw
from anoisetools.test.scripts import data_path

class Sff2RawTestCase(unittest.TestCase):

    def setUp(self):
        self.infile = data_path('test.sff')
        with open(data_path('test.raw')) as fp:
            self.expected = list(fp)[1:]

    def test_run(self):
        with tempfile.NamedTemporaryFile(prefix='test.raw', delete=False) as tf:
            try:
                sff2raw.sff_to_raw(self.infile, tf, 'test.raw')
                with open(tf.name) as tf:
                    v = list(tf)[1:]
            finally:
                os.remove(tf.name)
        self.assertEqual(self.expected, v)
