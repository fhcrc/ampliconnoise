
import cStringIO as StringIO
import unittest

from anoisetools.scripts import split

def _line_iterator(s):
    return iter(s.splitlines())

class ReadBarcodeTestCase(unittest.TestCase):

    def test_basic(self):
        s = "Barcode1,ACGTGA"
        i = _line_iterator(s)
        actual = split._load_barcodes(i)
        expected = {'ACGTGA': 'Barcode1'}
        self.assertEqual(expected, actual)

    def test_extra_data(self):
        """
        Extra items in a row should be ignored
        """
        s = """Barcode1,ACGTGA,Data from patient X
Barcode2,ACTGAA,Data from patient Y"""
        i = _line_iterator(s)
        actual = split._load_barcodes(i)
        expected = {'ACGTGA': 'Barcode1', 'ACTGAA':'Barcode2'}
        self.assertEqual(expected, actual)

    def test_no_dup_barcodes(self):
        """
        Duplicate barcodes are not allowed
        """
        s = "Barcode1,ACGTGA\nBarcode2,ACGTGA"
        i = _line_iterator(s)
        self.assertRaises(ValueError, split._load_barcodes, i)

    def test_no_dup_labels(self):
        """
        Duplicate labels are not allowed
        """
        s = "Barcode1,ACGTGA\nBarcode1,ACGTGG"
        i = _line_iterator(s)
        self.assertRaises(ValueError, split._load_barcodes, i)
