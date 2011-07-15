
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
        self.assertEquals([99, 2, 109, 7],
                          record0.annotations['flow_values'])

class AnoiseRecordConverterTestCase(unittest.TestCase):

    def setUp(self):
        flows = [104, 1, 102, 7, 5, 97, 6, 200, 96, 109, 98, 8, 2, 118]
        bases = sff.flow_to_seq(flows)
        sequence = Seq(bases, generic_dna)
        self.record = SeqRecord(sequence, id='FTWCYXX01BTPDQ')
        self.record.annotations['clip_qual_right'] = 14
        self.record.annotations['flow_values'] = flows

    def test_to_anoise_raw(self):
        record = self.record
        expected = ">FTWCYXX01BTPDQ\n14 1.04 0.01 1.02 0.07 0.05 0.97 0.06 2.00 0.96 1.09 0.98 0.08 0.02 1.18"
        self.assertEqual(expected, anoiseio._record_to_anoise_raw(record))
