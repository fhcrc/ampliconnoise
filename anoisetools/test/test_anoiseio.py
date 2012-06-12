
import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from anoisetools import anoiseio, sff


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
        self.assertEquals('GYUL5KN01A46DL', record0.id)
        self.assertEquals([99, 2, 109, 7, 6],
                          record0.annotations['flow_values'])
        self.assertEquals(4, record0.annotations['clip_flow_right'])

class AnoiseRecordConverterTestCase(unittest.TestCase):

    def setUp(self):
        flows = (101, 1, 99, 2, 0, 103, 2, 100, 0, 111, 110, 1, 3, 108, 0, 2,
                107, 0, 0, 103, 101, 1, 13, 102, 0, 8, 204, 12, 0, 105, 1, 290,
                2, 92, 7, 98, 9, 6, 98, 4, 101, 10, 0, 199, 5, 299, 3, 101, 4,
                99)
        bases = sff.flow_to_seq(flows)
        sequence = Seq(bases, generic_dna)
        self.record = SeqRecord(sequence, id='FTWCYXX01BTPDQ')
        self.record.annotations['clip_qual_right'] = 14
        self.record.annotations['clip_qual_left'] = 4
        self.record.annotations['flow_values'] = flows
        self.record.annotations['flow_index'] = (1, 2, 3, 2, 2, 1, 3, 3, 3, 1,
                3, 3, 0, 3, 2, 0, 0, 2, 2, 3, 2, 3, 0, 2, 0, 0, 2, 2, 3, 0, 2,
                1, 3, 3, 0, 3, 0, 1, 0, 3, 2, 0, 2, 3, 0, 3, 0, 2, 0, 0)

    def test_to_anoise_raw(self):
        record = self.record
        expected = ">FTWCYXX01BTPDQ\n30 1.01 0.01 0.99 0.02 0.00 1.03 0.02 1.00 0.00 1.11 1.10 0.01 0.03 1.08 0.00 0.02 1.07 0.00 0.00 1.03 1.01 0.01 0.13 1.02 0.00 0.08 2.04 0.12 0.00 1.05 0.01 2.90 0.02 0.92 0.07 0.98 0.09 0.06 0.98 0.04 1.01 0.10 0.00 1.99 0.05 2.99 0.03 1.01 0.04 0.99"
        self.assertEqual(expected, anoiseio._record_to_anoise_raw(record))
