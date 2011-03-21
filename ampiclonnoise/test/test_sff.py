
import unittest

from ampiclonnoise import sff

class TestFlowToSeq(unittest.TestCase):

    def setUp(self):
        self.test_flow = [0.18, 1.03, 1.02, 0.70, 1.12, 0.07, 0.14, 4.65]
        self.test_expected = 'ACGTGGGGG'

    def test_flow_to_seq(self):
        actual = sff.flow_to_seq(self.test_flow)
        self.assertEqual(self.test_expected, actual)

class TestFlowerReader(unittest.TestCase):

    def setUp(self):
        self.test_iter = """>FTWCYXX01BTPDQ
  Info: \t2009-04-08 11:31:35 R1 (631,1068)
  Clip: \t5 282
  Flows:\t1.04 0.01 1.02 0.07 0.05 0.97 0.06 2.00 0.96 1.09 0.98 0.08 0.02 1.18
  Index:\t1 3 6 8 8 9 10 11 14 16 18 18 18 18 18 21 21 23 23 23 23 25 27 27 27 3
  Bases:\ttcagGTACAGAAAAATTCCCCTCCCAATTAAAACTATGTGATGTGATTTCTATGTCCCCTCCTGAGGGTT
  Quals:\t37 37 37 37 37 37 37 37 38 34 20 20 20 20 15 34 34 37 37 37 37 38 34 3
>FTWCYXX01B0TTO
  Info: \t2009-04-08 11:31:35 R1 (712,1642)
  Clip: \t5 278
  Flows:\t1.04 0.04 1.01 0.09 0.06 0.99 0.09 1.97 1.02 1.06 0.98 0.07 0.01 1.20
  Index:\t1 3 6 8 8 9 10 11 14 16 18 18 18 18 18 21 21 23 25 27 27 29 31 33 35 3
  Bases:\ttcagGTACAGAAAAATTCTCCTCTCAATTAAAACTATGTGATGTGATTTCTATGTCCCCTCCTGAGGGTT
  Quals:\t37 37 37 37 37 37 37 37 38 38 28 28 28 28 28 38 38 37 37 37 37 37 40""".splitlines()

        self.reader = sff.parse_flower(iter(self.test_iter))

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

    def test_invalid_key(self):
        self.test_iter.insert(2, '  Invalid:\t0 0 0 0 0')
        reader = sff.parse_flower(iter(self.test_iter))
        self.assertRaisesRegexp(ValueError, 'Unknown key: Invalid',
                next, reader)

    def test_duplicate_key(self):
        # Duplicate a line
        self.test_iter.insert(2, self.test_iter[2])
        reader = sff.parse_flower(iter(self.test_iter))
        with self.assertRaisesRegexp(ValueError, '^Clip already set:') as cm:
            next(reader)

    def test_validate(self):
        self.test_iter.pop(2)
        self.test_iter.pop(3)
        reader = sff.parse_flower(iter(self.test_iter))
        with self.assertRaisesRegexp(ValueError, r'^Missing attribute\(s\):') as cm:
            next(reader)

class TestSFFRecord(unittest.TestCase):
    def setUp(self):
        self.text_record = """>FTWCYXX01BTPDQ
  Info: \t2009-04-08 11:31:35 R1 (631,1068)
  Clip: \t5 282
  Flows:\t1.04 0.01 1.02 0.07 0.05 0.97 0.06 2.00 0.96 1.09 0.98 0.08 0.02 1.18
  Index:\t1 3 6 8 8 9 10 11 14 16 18 18 18 18 18 21 21 23 23 23 23 25 27 27 27 3
  Bases:\ttcagGTACAGAAAAATTCCCCTCCCAATTAAAACTATGTGATGTGATTTCTATGTCCCCTCCTGAGGGTT
  Quals:\t37 37 37 37 37 37 37 37 38 34 20 20 20 20 15 34 34 37 37 37 37 38 34 3"""
        self.reader = sff.parse_flower(iter(self.text_record.splitlines()))
        self.record = next(self.reader)

    def test_str(self):
        actual = str(self.record)
        actual_lines = actual.splitlines()
        expected_lines = self.text_record.splitlines()

        self.assertEquals(len(expected_lines), len(actual_lines),
                          msg="Line counts should match")

        for expected_line, actual_line in zip(expected_lines, actual_lines):
            self.assertEquals(expected_line, actual_line)

    def test_to_anoise_raw(self):
        expected = ">FTWCYXX01BTPDQ\n14 1.04 0.01 1.02 0.07 0.05 0.97 0.06 2.00 0.96 1.09 0.98 0.08 0.02 1.18"
        self.assertEqual(expected, self.record.to_anoise_raw())
