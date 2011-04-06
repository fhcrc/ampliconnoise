import unittest

from anoisetools import sff

class TestFlowToSeq(unittest.TestCase):

    def setUp(self):
        self.test_flow = [0.18, 1.03, 1.02, 0.70, 1.12, 0.07, 0.14, 4.65]
        self.test_expected = 'ACGTGGGGG'
        self.test_flow_int = [int(i*100) for i in self.test_flow]

    def test_flow_to_seq(self):
        actual = sff.flow_to_seq(self.test_flow)
        self.assertEqual(self.test_expected, actual)

    def test_flow_to_seq_int(self):
        actual = sff.flow_to_seq(self.test_flow_int)
        self.assertEqual(self.test_expected, actual)

