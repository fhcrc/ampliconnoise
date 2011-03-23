"""
Given the input:

  * A flower-converted Standard Flowgram File
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

from ampliconnoise import sff, anoiseio


class _FlowerWriter(object):
    def __init__(self, fp, *args):
        self._fp = fp

    def write(self, record):
        print >>self._fp, str(record)

    def close(self):
        self._fp.close()

# Output formatters - each takes a record, returns a string to write
WRITERS = {'flower': _FlowerWriter, 'anoise_raw': anoiseio.AnoiseRawWriter }


def build_parser(subparsers):
    """
    Adds description, command line options to parser
    """
    parser = subparsers.add_parser('split', help="Split .sff.txt by tag")
    parser.epilog = __doc__
    parser.description = "sff.txt splitter"
    parser.formatter_class = argparse.RawDescriptionHelpFormatter

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
    parser.add_argument('--output-directory', metavar='DIR', default='.',
            help='Output directory for split files (default: %(default)s)')
    parser.add_argument('--output-format', metavar='FORMAT',
            default='anoise_raw', choices=WRITERS.keys(),
            help="Output format (choices: [%(choices)s], default: %(default)s)")
    parser.add_argument('--unmatched-name', default='unmatched',
            help='Name for file containing unmatched records. '
                 'Default: %(default)s)')

    return parser

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
    """
    Splits flows based on the initial sequence contained in the flowgram
    """

    def __init__(self, barcode_map, primer, dest_dir, unmatched_dest,
                 writer_cls):
        """
        Initialize a new instance

        :param map barcode_map: map from sequence barcode to name
        :param str primer: Primer sequence
        :param str dest_dir: destination directory for output
        :param str unmatched_dest: name of file for all input
          sequences not matching any barcodes.
        :param func writer: Function to use to format records. Passed a
          single argument as input - the record
        """
        self.barcode_map = barcode_map
        self.primer = primer
        self.dest_dir = dest_dir
        self.unmatched_dest = unmatched_dest
        self.writer_cls = writer_cls

        # Compile barcode matching Regex
        min_barcode_length = min(len(k) for k in self.barcode_map)
        max_barcode_length = max(len(k) for k in self.barcode_map)
        self._barcode_re = re.compile(r'^{0}(\w{{{1},{2}}}){3}'.format(
            sff.READ_BEGIN, min_barcode_length, max_barcode_length, primer),
                                      re.IGNORECASE)
        _makedirs(dest_dir)
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
            fp = open(outpath, 'w')
            self._handles[barcode] = self.writer_cls(fp, barcode)

        # Create a default handle
        default_outpath = os.path.join(self.dest_dir, self.unmatched_dest)
        fp = open(default_outpath + '.raw', 'w')
        self._handles[None] = self.writer_cls(fp, 'Unknown')

    def close(self):
        """
        Closes all open file handles on this instance
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

        writer = self._handles.get(barcode, self._handles[None])
        writer.write(record)

        return barcode

    def split(self, iterable):
        """
        Takes an iterable generating :ref:`ampliconnoise.sff.SFFRead`s,
        writes them to a set of output files
        returns a dictionary of mapping barcode -> # of reads
        """
        counts = collections.defaultdict(int)
        unmatched_counts = collections.defaultdict(int)

        for record in iterable:
            barcode = self._handle_record(record)
            counts[barcode] += 1

        if barcode in self._handles:
            counts[barcode] += 1
        else:
            unmatched_counts[barcode] += 1

        return counts


def main(parsed_args):
    """
    Command-line functionality

    * parse arguments
    * read barcodes
    * run
    """
    writer = WRITERS[parsed_args.output_format]

    # Split barcodes
    with parsed_args.barcode_file:
        barcodes = _load_barcodes(parsed_args.barcode_file)

    splitter = SFFRunSplitter(barcodes, parsed_args.primer, parsed_args.output_directory,
                              parsed_args.unmatched_name, writer)

    # Run
    with contextlib.closing(splitter):
        splitter.open()
        with parsed_args.sff_file:
            reader = sff.parse_flower(parsed_args.sff_file)
            result = splitter.split(reader)

    # Special treatment for unmatched records
    unmatched = result[None]
    del result[None]

    items = sorted(result.items())
    for k, v in items:
        print 'Barcode {0}: {1:5d} records'.format(k, v)

    print 'Unmatched:', unmatched, 'records'
