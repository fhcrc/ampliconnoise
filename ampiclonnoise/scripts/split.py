"""
``split``
=========

.sff.txt splitter by barcode.

Given the input:

  * A text-converted Standard Flowgram File
  * A file containing mappings from DNA sequence barcodes to a meaningful
    identifier
  * A regular expression to identify the primer sequence used

Writes to a series of files named after the meaningful identifiers provided,
each containing only those sequences from the input sff with the matching
barcode.

Any unmatched sequences are written to a separate file
"""

import argparse
import collections
import contextlib
import csv
import errno
import os
import os.path
import re
import sys

from ampiclonnoise import sff


def _makedirs(d):
    """
    Make all directories to leaf dir ``d``
    """
    try:
        os.makedirs(d)
    except OSError, e:
        if not e.errno == errno.EEXIST:
            raise e


def _load_barcodes(fp):
    """
    Read barcodes from fp

    Raises a ValueError on duplicate key *or* value.

    fp should contain (at least) sequence_identifier,barcode_sequence
    pairs.
    """
    reader = csv.reader(fp)
    values = set()
    d = {}
    for line in reader:
        k, v = line[1], line[0]
        if k in d:
            raise ValueError("Duplicate barcode: {0}".format(k))
        if v in values:
            raise ValueError("Duplicate value: {0}".format(v))
        d[k] = v
        values.add(v)
    return d

def _close_all(files):
    """
    Closes an iterable of files, aggregating any exceptions
    """
    exceptions = []
    for f in files:
        try:
            f.close()
        except Exception, e:
            exceptions.append(e)

    if exceptions:
        raise IOError("Could not close {0} files: {1}".format(
            len(exceptions), exceptions))


class SFFRunSplitter(object):

    def __init__(self, barcode_map, primer, dest_dir, unmatched_dest):
        """
        Initialize a new instance

        :param map barcode_map: map from sequence barcode to name
        :param str primer: Primer sequence
        :param str dest_dir: destination directory for output
        :param str unmatched_dest: name of file for all input
          sequences not matching any barcodes.
        """
        self.barcode_map = barcode_map
        self.primer = primer
        min_barcode_length = min(len(k) for k in self.barcode_map)
        max_barcode_length = max(len(k) for k in self.barcode_map)
        self._barcode_re = re.compile(r'TCAG(\w{{{0},{1}}}){2}'.format(
            min_barcode_length, max_barcode_length, primer),
                                      re.IGNORECASE)
        self.dest_dir = dest_dir
        _makedirs(dest_dir)
        self.unmatched_dest = unmatched_dest
        self._handles = None

    def open(self):
        """
        Opens file handles for each barcode, in preparation for writing.
        """
        if self._handles:
            raise IOError("Already initialized with {0} handles".format(
                     len(self._handles)))
        self._handles = {}

        for barcode, name in self.barcode_map.items():
            fname = name + '.raw'
            outpath = os.path.join(self.dest_dir, fname)
            self._handles[barcode] = open(outpath, 'w')

        # Create a default handle
        default_outpath = os.path.join(self.dest_dir, self.unmatched_dest)
        self._handles[None] = open(default_outpath + '.raw', 'w')

    def close(self):
        """
        Closes all open file handles
        """
        if self._handles is not None:
            _close_all(self._handles.values())

    def _handle_record(self, record):
        """
        Identifies the barcode in the record, writes the record
        to the appropriate outfile, and returns the barcode
        """
        bases = record.bases_from_flows()
        m = self._barcode_re.match(bases)
        barcode = m.group(1) if m else None
        fp = self._handles.get(barcode, self._handles[None])
        print >> fp, record
        return barcode

    def split(self, iterable):
        """
        Takes an iterable generating :ref:`ampiclonnoise.sff.SFFRead`s,
        writes them to a set of output files
        returns a dictionary of mapping barcode -> # of reads
        """
        counts = collections.defaultdict(int)
        for record in iterable:
            barcode = self._handle_record(record)
            counts[barcode] += 1

        return counts


def main(args=sys.argv[1:]):
    """
    Command-line functionality

    * parse arguments
    * read barcodes
    * run
    """
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument('barcode_file', metavar='BARCODE_FILE',
            type=argparse.FileType('r'),
            help='Path to barcode file with barcode_id,barcode_seq pairs')
    parser.add_argument('primer', metavar='PRIMER',
            help='Regular expression identifying the primer used')

    # Optional arguments
    parser.add_argument('--sff-file', metavar='SFF_TXT', default=sys.stdin,
            type=argparse.FileType('r'),
            help='Flower-decoded SFF text file. Default: stdin')
    parser.add_argument('--output-directory', metavar='DIR',
            default='split', help='Output directory for split files')
    parser.add_argument('--unmatched-name', default='unmatched',
            help='Name for file containing unmatched records. '
                 'Default: %(default)s)')

    parsed = parser.parse_args(args)

    # Split barcodes
    with parsed.barcode_file:
        barcodes = _load_barcodes(parsed.barcode_file)

    splitter = SFFRunSplitter(barcodes, parsed.primer, parsed.output_directory,
                              parsed.unmatched_name)
    with contextlib.closing(splitter):
        splitter.open()
        with parsed.sff_file:
            reader = sff.parse_flower(parsed.sff_file)
            result = splitter.split(reader)

    # Special treatment for unmatched records
    unmatched = result[None]
    del result[None]

    items = sorted(result.items())
    for k, v in items:
        print 'Barcode {0}: {1:5d} records'.format(k, v)
    print 'Unmatched:', unmatched, 'records'