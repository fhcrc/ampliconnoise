import os
import os.path
import shutil
import unittest

from anoisetools.scripts import dispatch
from anoisetools.test.scripts import LineCompareMixIn, data_path, tempdir, cd

is_enabled = os.environ.get('TEST_PYRONOISE') is not None

@unittest.skipUnless(is_enabled, 'Set environment variable TEST_PYRONOISE')
class PyroNoiseTestCase(LineCompareMixIn, unittest.TestCase):
    def setUp(self):
        self.sff_path = data_path('test.sff')
        self.expected_path = data_path('test-pnoise_cd.fa',)
        self.mapping_path = data_path('test-pnoise.mapping')
        with open(self.expected_path) as fp:
            self.expected = list(fp)

    def test_run_pyronoise(self):
        sff_path = os.path.abspath(self.sff_path)
        sff_name = os.path.basename(sff_path)
        args = ['pyronoise', sff_name, '--mpi-args', '-np 2', '--stub', 'test']
        with tempdir(prefix='pyronoise-') as td:
            with cd(td):
                self.assertEqual(td, os.getcwd())
                shutil.copy(sff_path, sff_name)
                dispatch.main(args)
                self.assertEqual(td, os.getcwd())
            try:
                self.assertLinesEqual(self.expected_path, 'test-pnoise_cd.fa')
                self.assertLinesEqual(self.mapping_path, 'test-pnoise.mapping')
            finally:
                for f in ('test-pnoise_cd.fa', 'test-pnoise.mapping'):
                    if os.path.exists(f):
                        os.remove(f)
