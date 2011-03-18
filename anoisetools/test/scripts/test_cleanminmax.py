"""
Tests for cleanminmax.py
"""

from cStringIO import StringIO
import os.path
import re
import unittest

from anoisetools import flower
from anoisetools.scripts import cleanminmax


class TestFlowgramValid(unittest.TestCase):
    def test_basic(self):
        flowgram = [1.0, 1.0, 0.0, 1.0]
        self.assertEquals(True, cleanminmax.is_flowgram_valid(flowgram))

    def test_high_signal(self):
        flowgram_accept = [1.0, 1.0, 9.49, 0.2]
        flowgram_reject = [1.0, 1.0, 9.51, 0.2]
        self.assertEquals(True, cleanminmax.is_flowgram_valid(flowgram_accept))
        self.assertEquals(False, cleanminmax.is_flowgram_valid(flowgram_reject))

    def test_low_signal(self):
        flowgram_accept = [1.0, 1.0, 1.0, 0.7]
        flowgram_reject = [1.0, 1.0, 1.0, 0.6]
        self.assertEquals(True, cleanminmax.is_flowgram_valid(flowgram_accept))
        self.assertEquals(False, cleanminmax.is_flowgram_valid(flowgram_reject))

    def test_no_signal(self):
        flowgram_nosignal = [0.0] * 4
        self.assertEquals(False, cleanminmax.is_flowgram_valid(flowgram_nosignal))

    def test_arg_count(self):
        self.assertRaisesRegexp(ValueError, "^Unexpected flowgram length",
                cleanminmax.is_flowgram_valid, range(5))

class TestHandleRecord(unittest.TestCase):

    def setUp(self):
        self.flows = [1.03, 0.00, 1.05, 0.01, 0.03, 0.99, 0.05, 0.93, 0.04,
                0.95, 0.04, 0.99, 1.00, 0.00, 0.99, 0.03, 1.01, 0.00, 0.96,
                0.01, 0.06, 2.93, 0.04, 0.01, 2.08, 0.03, 0.05, 2.93, 0.09,
                0.00, 2.05, 0.09, 1.04, 0.04, 0.08, 1.03, 0.05, 4.00, 0.07,
                0.03, 0.98, 0.11, 2.02, 0.01, 0.08, 1.06, 0.10, 0.00, 1.00,
                1.03, 0.98, 0.01, 0.07, 1.99, 0.10, 0.00, 1.02, 1.00, 1.02,
                0.00, 0.96, 0.09, 1.96, 0.04, 0.10, 0.99, 0.15, 0.95, 1.00,
                0.07, 0.05, 1.02, 2.98, 0.00, 0.06, 0.98, 0.10, 0.02, 1.09,
                0.05, 0.97, 1.11, 0.09, 0.02, 1.01, 3.13, 0.05, 1.01, 0.08,
                4.90, 0.05, 1.87, 0.06, 0.94, 0.06, 0.08, 0.98, 0.97, 0.07,
                0.97, 0.12, 0.13, 1.03, 0.09, 0.05, 0.98, 1.03, 0.03, 1.00,
                3.03, 0.14, 0.08, 1.00, 0.12, 0.05, 2.09, 0.12, 1.02, 0.04,
                0.87, 0.08, 3.98, 0.07, 0.15, 2.00, 1.03, 0.06, 1.01, 0.98,
                0.95, 0.05, 1.00, 0.07, 0.94, 0.05, 0.16, 4.21, 0.95, 0.07,
                0.98, 0.08, 1.01, 0.05, 1.03, 0.09, 0.93, 0.08, 1.05, 0.06,
                0.11, 1.02, 0.08, 0.96, 0.11, 0.97, 0.05, 0.11, 2.05, 0.13,
                0.06, 0.95, 3.74, 0.11, 0.97, 0.12, 2.07, 1.00, 0.08, 0.09,
                1.02, 0.93, 0.11, 0.05, 0.90, 0.12, 2.03, 0.09, 0.97, 0.99,
                0.15, 2.06, 0.10, 0.98, 0.15, 1.07, 0.14, 0.10, 3.10, 0.13,
                1.92, 0.13, 0.83, 2.01, 0.11, 1.02, 0.11, 0.13, 2.10, 0.18,
                0.13, 2.01, 0.98, 0.05, 1.95, 0.14, 1.94, 0.08, 0.10, 1.04,
                1.07, 4.09, 0.11, 0.11, 0.95, 0.16, 0.05, 0.94, 0.11, 2.02,
                0.08, 0.96, 0.11, 0.16, 0.97, 0.08, 0.11, 1.01, 0.07, 0.05,
                1.05, 0.14, 1.98, 2.94, 5.40, 0.13, 0.83, 0.13, 5.57, 0.15,
                0.10, 0.91, 0.25, 0.98, 0.12, 0.12, 1.13, 0.13, 1.02, 0.95,
                1.98, 1.01, 0.12, 0.13, 1.02, 0.12, 1.05, 0.98, 0.99, 1.01,
                0.10, 0.87, 1.03, 0.09, 1.07, 0.13, 1.02, 0.10, 0.07, 0.95,
                0.15, 0.09, 1.08, 1.04, 0.09, 0.08, 4.82, 0.14, 0.92, 0.08,
                0.16, 1.04, 0.12, 0.14, 1.10, 0.09, 0.07, 2.01, 0.11, 0.94,
                1.08, 0.14, 0.10, 4.99, 0.11, 1.03, 0.12, 0.09, 0.95, 0.15,
                0.88, 2.11, 0.12, 2.07, 0.15, 2.94, 1.01, 0.07, 0.99, 0.08,
                0.85, 0.14, 0.07, 1.15, 2.98, 0.15, 1.04, 0.04, 0.90, 0.12,
                0.10, 3.86, 0.90, 0.07, 0.98, 0.16, 3.08, 0.09, 0.95, 0.98,
                0.80, 0.12, 0.10, 1.05, 1.13, 1.05, 0.06, 0.94, 0.16, 0.10,
                1.00, 0.17, 0.16, 1.00, 0.06, 0.06, 1.16, 0.14, 0.09, 2.05,
                0.11, 1.02, 0.06, 0.09, 0.98, 1.93, 0.05, 0.09, 0.97, 0.06,
                0.04, 1.02, 0.94, 2.09, 0.09, 0.89, 1.01, 0.08, 1.05, 1.03,
                1.00, 0.10, 0.00, 0.93, 3.00, 1.04, 0.11, 0.09, 2.00, 0.06,
                0.10, 1.13, 0.15, 0.16, 1.12, 0.06, 1.04, 0.15, 0.93, 0.12,
                1.16, 1.21, 0.10, 0.09, 0.99, 1.98, 0.09, 0.10, 1.05, 0.08,
                1.98, 0.10, 1.89, 2.03, 0.14, 0.09, 0.99, 0.14, 2.01, 0.88,
                1.20, 0.14, 0.12, 1.03, 0.16, 1.08, 0.05, 0.14, 1.01, 0.09,
                1.08, 1.11, 1.05, 0.97, 0.09, 0.14, 2.12, 0.15, 0.10, 0.88,
                0.21, 0.13, 0.99, 0.98, 0.14, 0.17, 0.97, 0.15, 0.15, 0.93,
                0.13, 2.04, 0.11, 1.92, 0.09, 0.17, 1.06, 1.06, 0.06, 0.11,
                1.97, 0.16, 1.95, 0.09, 1.10, 0.12, 0.08, 1.08, 0.10, 0.19,
                2.08, 0.14, 3.06, 0.11, 2.03, 0.08, 1.06, 0.12, 0.10, 1.01,
                0.16, 0.88, 0.22, 0.15, 1.07, 1.97, 0.11, 0.11, 1.13, 0.13,
                0.95, 0.16, 0.13, 1.06, 0.08, 0.09, 2.02, 0.15, 0.08, 0.93,
                1.20, 0.11, 0.15, 2.01, 0.15, 2.13, 0.04, 0.14, 1.03, 0.07,
                0.89, 0.08, 1.01, 0.20, 0.94, 1.09, 0.96, 0.14]

    def tearDown(self):
        pass

    def test_min_length(self):
        assert len(self.flows) >= 400

        regex = re.compile('^TCAG.*?(A.*)')
        actual = cleanminmax.handle_record(self.flows, regex, 400, 600)
        self.assertTrue(actual[0] is not None)
        self.assertTrue(actual[1] is not None)

        self.fail("Check results")

class PerlOutputCompare(unittest.TestCase):

    def setUp(self):
        d = os.path.dirname(__file__)
        infile = os.path.join(d, 'data', 'AA-1-pol756.raw')
        self.fp = open(infile)
        # Ignoring first line
        self.fp.readline()
        with open(infile + '.fa') as fp:
            self.fasta_compare = fp.read()
        with open(infile + '.dat') as fp:
            self.dat_compare = fp.read()

    def tearDown(self):
        self.fp.close()

    def test_invoke(self):
        def reader(fp):
            while True:
                header = next(fp)
                header = header[1:]
                data = next(fp)
                data = map(float, data.split()[1:])
                record = flower.FlowerRecord(header)
                record.flows = data
                yield record
        fasta_handle = StringIO()
        dat_handle = StringIO()
        primer = 'AATTGGGCCTGAAAATCCATA'
        min_length = cleanminmax.DEFAULT_MIN_LENGTH
        max_length = 720

        cleanminmax.invoke(reader(self.fp), fasta_handle, dat_handle,
                           primer, min_length, max_length)

        fasta_str = fasta_handle.getvalue()
        dat_str = dat_handle.getvalue()

        self.assertTrue(len(fasta_str) > 0)
        self.assertTrue(len(dat_str) > 0)

        self.assertEquals(self.fasta_compare.splitlines(), fasta_str.splitlines())
        self.assertEquals(self.dat_compare, dat_str)

        self.fail("todo")


