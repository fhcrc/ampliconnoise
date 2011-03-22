
import unittest

from ampliconnoise import anoiseio, sff


class SplitKeysReaderTestCase(unittest.TestCase):

    def setUp(self):
        raw_string = '''0 CAGTCGT AA-1-pol756.raw
>GYUL5KN01A46DL
4 0.99 0.02 1.09 0.07 0.06
>GYUL5KN01EOUKP
12 1.04 0.06 1.05 0.07 0.09 0.99 0.09 1.08 0.08 0.02 1.05 0.03 0.02'''
        self.fp = iter(raw_string.splitlines())
        self.reader = anoiseio.AnoiseRawReader(self.fp)

    def test_header(self):
        reader = self.reader
        self.assertEquals(0, reader.count)
        self.assertEquals('CAGTCGT', reader.barcode)
        self.assertEquals('AA-1-pol756.raw', reader.fname)

    def test_record_count(self):
        records = [i for i in self.reader]
        self.assertEquals(2, len(records))

    def test_record(self):
        records = [i for i in self.reader]

        record0 = records[0]
        self.assertEquals('GYUL5KN01A46DL', record0.identifier)
        self.assertEquals([0.99, .02, 1.09, .07], record0.flows)

class AnoiseRecordConverterTestCase(unittest.TestCase):

    def setUp(self):
        self.record = sff.SFFRead('FTWCYXX01BTPDQ')
        self.record.clip = [0, 14]
        self.record.flows = [1.04, 0.01, 1.02, 0.07, 0.05, 0.97, 0.06, 2.00, 0.96, 1.09, 0.98, 0.08, 0.02, 1.18]

    def test_to_anoise_raw(self):
        record = self.record
        expected = ">FTWCYXX01BTPDQ\n14 1.04 0.01 1.02 0.07 0.05 0.97 0.06 2.00 0.96 1.09 0.98 0.08 0.02 1.18"
        self.assertEqual(expected, anoiseio._record_to_anoise_raw(record))
