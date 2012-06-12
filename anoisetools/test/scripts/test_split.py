import unittest

from anoisetools.scripts import split

def _line_iterator(s):
    return iter(s.splitlines())

class ReadBarcodeTestCase(unittest.TestCase):

    def test_basic(self):
        s = "Barcode1,ACGTGA,ACGT"
        i = _line_iterator(s)
        actual = split.load_barcodes(i)
        self.assertEqual(1, len(actual))
        self.assertEqual(1, len(actual['ACGTGA']))
        self.assertEqual('Barcode1', actual['ACGTGA'][0][1])

    def test_extra_data(self):
        """
        Extra items in a row should be ignored
        """
        s = """Barcode1,ACGTGA,ACGT,Data from patient X
Barcode2,ACTGAA,ACGT,Data from patient Y"""
        i = _line_iterator(s)
        actual = split.load_barcodes(i)
        self.assertEqual(2, len(actual))
        self.assertEqual('Barcode1', actual['ACGTGA'][0][1])

    def test_no_dup_barcodes(self):
        """
        Duplicate barcodes are not allowed
        """
        s = "Barcode1,ACGTGA,ACGT\nBarcode2,ACGTGA,ACGT"
        i = _line_iterator(s)
        self.assertRaises(ValueError, split.load_barcodes, i)

    def test_no_dup_labels(self):
        """
        Duplicate labels are not allowed
        """
        s = "Barcode1,ACGTGA,ACGT\nBarcode1,ACGTGG,ACGT"
        i = _line_iterator(s)
        self.assertRaises(ValueError, split.load_barcodes, i)
