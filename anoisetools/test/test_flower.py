
import unittest

from anoisetools import flower

class TestFlowToSeq(unittest.TestCase):

    def setUp(self):
        self.test_flow = [0.18, 1.03, 1.02, 0.70, 1.12, 0.07, 0.14, 4.65]
        self.test_expected = 'ACGTGGGGG'

    def test_flow_to_seq(self):
        actual = flower.flow_to_seq(self.test_flow)
        self.assertEqual(self.test_expected, actual)

class TestFlowerReader(unittest.TestCase):

    def setUp(self):
        self.test_iter = """>FTWCYXX01BTPDQ
  Info:   2009-04-08 11:31:35 R1 (631,1068)
  Clip:   5 282
  Flows:  1.04 0.01 1.02 0.07 0.05 0.97 0.06 2.00 0.96 1.09 0.98 0.08 0.02 1.18
  Index:  1 3 6 8 8 9 10 11 14 16 18 18 18 18 18 21 21 23 23 23 23 25 27 27 27 3
  Bases:  tcagGTACAGAAAAATTCCCCTCCCAATTAAAACTATGTGATGTGATTTCTATGTCCCCTCCTGAGGGTT
  Quals:  37 37 37 37 37 37 37 37 38 34 20 20 20 20 15 34 34 37 37 37 37 38 34 3
  >FTWCYXX01B0TTO
  Info:   2009-04-08 11:31:35 R1 (712,1642)
  Clip:   5 278
  Flows:  1.04 0.04 1.01 0.09 0.06 0.99 0.09 1.97 1.02 1.06 0.98 0.07 0.01 1.20
  Index:  1 3 6 8 8 9 10 11 14 16 18 18 18 18 18 21 21 23 25 27 27 29 31 33 35 3
  Bases:  tcagGTACAGAAAAATTCTCCTCTCAATTAAAACTATGTGATGTGATTTCTATGTCCCCTCCTGAGGGTT
  Quals:  37 37 37 37 37 37 37 37 38 38 28 28 28 28 28 38 38 37 37 37 37 37 40""".splitlines()

        self.reader = flower.read_flower(iter(self.test_iter))

    def tearDown(self):
        pass

    def test_read_1(self):
        record = self.reader.next()
        self.assertEquals("2009-04-08 11:31:35 R1 (631,1068)", record.info)
        self.assertEquals([5, 282], record.clip)

        self.assertEquals(14, len(record.flows))

    def test_read_clip(self):
        record = self.reader.next()
        self.assertEquals([5, 282], record.clip)
        record = self.reader.next()
        self.assertEquals([5, 278], record.clip)

    def test_read_flows(self):
        record = self.reader.next()
        self.assertEquals(14, len(record.flows))
        self.assertEquals(1.04, record.flows[0])
        self.assertEquals(1.18, record.flows[13])

    def test_read_identifier(self):
        record = self.reader.next()
        self.assertEquals("FTWCYXX01BTPDQ", record.identifier)

    def test_read_bases(self):
        record = self.reader.next()
        self.assertEquals('tcagGTACAGAAAAATTCCCCTCCCAATTAAAACTATGTGATGTGATTTCTATGTCCCCTCCTGAGGGTT',
                          record.bases)
        record = self.reader.next()
        self.assertEquals('tcagGTACAGAAAAATTCTCCTCTCAATTAAAACTATGTGATGTGATTTCTATGTCCCCTCCTGAGGGTT',
                          record.bases)

    def test_read_all(self):
        """
        Make sure all the records are read
        """
        records = list(self.reader)
        self.assertEquals(2, len(records))
        self.assertEquals(['FTWCYXX01BTPDQ', 'FTWCYXX01B0TTO'],
                          [record.identifier for record in records])
