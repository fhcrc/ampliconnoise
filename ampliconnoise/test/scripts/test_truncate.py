from cStringIO import StringIO
import unittest

from ampliconnoise.scripts import truncate

class TrimTestCase(unittest.TestCase):

    def setUp(self):
        self.input = iter('''>GYUL5KN02IY0NH qclip: 5..519
GCGATCTAAGTTCATAACCCATCCAAGGAATGGAGGTTCTTTCTGATGCTTTTTGTCTGGTGTAGTGAATCCCCAGCTCAATAGGTGAGCT'''.splitlines())
        self.barcode = 'GCGATCT'
        self.output = StringIO()

    @property
    def _value(self):
        return self.output.getvalue()

    def test_barcode_removed(self):
        truncate.trim(self.barcode, 20, self.input, self.output)
        self.assertFalse(self._value.splitlines()[1].startswith(self.barcode))

    def test_expected_full(self):
        truncate.trim(self.barcode, 20, self.input, self.output)
        self.assertEquals('''>GYUL5KN02IY0NH qclip: 5..519
AAGTTCATAACCCATCCAAG\n''', self._value)

