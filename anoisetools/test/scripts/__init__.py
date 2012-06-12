import contextlib
import os
import os.path
import shutil
import tempfile

def data_path(*args):
    return os.path.join(os.path.dirname(__file__), 'data', *args)

class LineCompareMixIn(object):
    def assertLinesEqual(self, path1, path2):
        def lines(path):
            with open(path) as fp:
                return list(i.rstrip() for i in fp)
        expected_lines = lines(path1)
        actual_lines = lines(path2)
        self.assertEqual(len(expected_lines), len(actual_lines))
        for e, a in zip(expected_lines, actual_lines):
            self.assertEqual(e, a)

@contextlib.contextmanager
def tempdir(**kwargs):
    td = tempfile.mkdtemp(**kwargs)
    try:
        yield td
    finally:
        shutil.rmtree(td)

@contextlib.contextmanager
def cd(path):
    orig_path = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(orig_path)
