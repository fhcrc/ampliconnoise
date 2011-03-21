import os.path
import shutil
import tempfile
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


class SFFRunSplitterTestCase(unittest.TestCase):

    def setUp(self):
        # Make an output dir for the output
        self.output_dir = tempfile.mkdtemp(prefix='tmp_anoisetest')
        self.barcodes = {'AAAA': 'A', 'CCCC': 'B'}
        self.primer = 'GGGGG[GC]GG'
        self.unmatched_dest = 'unmatched'
        self.instance = split.SFFRunSplitter(self.barcodes, self.primer,
                                             self.output_dir,
                                             self.unmatched_dest)

    def tearDown(self):
        self.instance.close()
        shutil.rmtree(self.output_dir)

    def test_barcode_re(self):
        instance = self.instance
        barcode_re = instance._barcode_re

        test = 'TCAGAAAAGGGGGGGGACGC'
        m = barcode_re.match(test)
        self.assertTrue(m is not None, '{0} must match'.format(test))
        self.assertEquals('AAAA', m.group(1),
                          msg='Barcode group must be captured')

        test = 'TGACAAAAGGGCGGGG'
        self.assertTrue(barcode_re.match(test) is None,
                        msg='Invalid initial sequence must not match')

    def test_open_files_present(self):
        instance = self.instance
        def check_name(key, fname):
            self.assertTrue(key in instance._handles, msg='Missing key: {0}'.format(key))
            path = instance._handles[key].name
            self.assertEquals(fname, os.path.basename(path))
        instance.open()
        self.assertEquals(3, len(instance._handles))
        check_name(None, 'unmatched.raw')
        check_name('AAAA', 'A.raw')
        check_name('CCCC', 'B.raw')

    def test_open_files_open(self):
        # Make sure all files are open
        instance = self.instance
        instance.open()
        for f in instance._handles.values():
            self.assertTrue(not f.closed, msg='Closed: {0}'.format(f))

    def test_open_files_close(self):
        instance = self.instance
        instance.open()
        instance.close()
        for f in instance._handles.values():
            self.assertTrue(f.closed, msg='Still open: {0}'.format(f))
