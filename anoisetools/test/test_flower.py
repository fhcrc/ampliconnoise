
import unittest

from anoisetools import flower

class TestFlowToSeq(unittest.TestCase):

    def setUp(self):
        self.test_flow = [0.18, 1.03, 1.02, 0.70, 1.12, 0.07, 0.14, 4.65]
        self.test_expected = 'ACGTGGGGG'

    def test_flow_to_seq(self):
        actual = flower.flow_to_seq(self.test_flow)
        self.assertEqual(self.test_expected, actual)

