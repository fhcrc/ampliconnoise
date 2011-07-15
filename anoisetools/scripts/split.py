"""
Split input FASTA files
"""

import argparse
import collections
import csv
import logging
import multiprocessing
import os
import os.path
import Queue
import re

from Bio import SeqIO


class SequenceWriter(object):
    """
    SequenceWriter writes sequences.

    Designed to be run as a subprocess, instances are initialized with a
    queue from which to draw sequences, an output file path, and an output
    format.

    Records are read from the queue and written to the output file until
    None is encountered in the queue.
    """

    def __init__(self, queue, path, sequence_format):
        self.count = multiprocessing.Value('i', 0)
        self.queue = queue
        self.path = path
        self.sequence_format = sequence_format
        self.cancel = False

    def __iter__(self):
        """
        Yield sequences from the input queue
        """
        while True:
            if self.cancel:
                break
            try:
                # Try to pull a new item, waiting two seconds
                record = self.queue.get(True, 2)
                if record is None:
                    break
                self.count.value += 1
                yield record
            except Queue.Empty:
                pass

    def write(self):
        """
        Write sequences passed to the input queue to self.path
        """
        SeqIO.write(iter(self), self.path, self.sequence_format)


def load_barcodes(fp):
    d = {}
    reader = csv.reader(fp)
    for record in reader:
        value, key = record[:2]
        if key in d:
            raise ValueError("Duplicate key: {0}".format(key))
        if value in d.values():
            raise ValueError("Duplicate value: {0}".format(value))
        d[key] = value

    return d

class AlreadyOpenError(BaseException):
    """
    Indicates that the object was already opened
    """
    pass


class SequenceSplitter(object):
    """
    Splits sequence files by barcode
    """
    barcode_start = 4

    def __init__(self, primer, barcodes, make_dirs=False,
            unmatched_name="unknown.sff", output_format='sff'):
        self.barcodes = barcodes

        min_len = min(len(barcode) for barcode in barcodes)
        max_len = max(len(barcode) for barcode in barcodes)
        self.primer = primer
        self.barcode_re = re.compile(
                '^(.{{{0},{1}}}){2}'.format(min_len, max_len, primer),
                re.IGNORECASE)
        self._output_files = {}
        self.output_format = output_format
        self.make_dirs = make_dirs
        self.unmatched_name = unmatched_name

    def open(self):
        """
        Create a subprocess running a SequenceWriter for each barcode.
        """
        if self._output_files:
            raise AlreadyOpenError("SequenceSplitter has already been opened.")

        writer_templates = []

        # Build a list of barcode, filepath tuples to construct writers from
        for barcode, identifier in self.barcodes.items():
            path = identifier + '.' + self.output_format
            if self.make_dirs:
                os.makedirs(identifier)
                path = os.path.join(identifier, path)

            writer_templates.append((barcode, path))

        # Special case: unmatched barcodes
        writer_templates.append((None, self.unmatched_name))

        for barcode, path in writer_templates:
            queue = multiprocessing.Queue()
            writer = SequenceWriter(queue, path, self.output_format)

            # Create a process to write output records
            process = multiprocessing.Process(target=writer.write)
            process.start()
            self._output_files[barcode] = queue, writer, process

    def close(self):
        """
        Close all the output files
        """
        for queue, writer, process in self._output_files.values():
            queue.put(None)
            process.join()

        self._output_files = {}

    def run(self, sequences):
        """
        Iterates over sequences, writing them to the file associated with their
        barcode.

        Any sequences which cannot be matched are written to a file named after
        unmatched_name
        """
        counter = collections.Counter()
        for sequence in sequences:
            seq_str = str(sequence.seq)[self.barcode_start:]
            m = self.barcode_re.match(seq_str)
            if m:
                barcode = m.group(1)
                counter[barcode] += 1
                queue = self._output_files.get(barcode,
                                               self._output_files[None])[0]
            else:
                queue = self._output_files[None][0]

            # There's a problem pickling BioPython _RestrictedDict
            # instances. Here's a hack to try to work around:
            # Convert the internal SeqRecord._per_letter_annotations
            # to a plain old Python dict.
            sequence._per_letter_annotations = sequence.letter_annotations.copy()

            queue.put(sequence)

        logging.info("%d sequences", sum(counter.values()))
        logging.info("Most common: %s" % '\n'.join(
                     ': '.join(map(str, i)) for i in counter.most_common(30)))


def build_parser(subparsers):
    """
    Adds description, command line options to parser
    """
    logging.basicConfig(level=logging.INFO)
    parser = subparsers.add_parser('split', help="Split an .sff")

    parser.add_argument("primer", help="""Sequence primer used,
            as regular expression""")
    parser.add_argument("barcode_file", type=argparse.FileType('r'),
            help="""Path to barcode file with barcode_id,barcode_seq pairs""")
    parser.add_argument("input_file", metavar="SFF",
            type=argparse.FileType('rb'),
            help="""Path to input SFF (use - for stdin)""")
    parser.add_argument("--make-dirs", action="store_true", default=False,
            help="Make a subdirectory for each output file")
    parser.add_argument('--unmatched-name',
            help="""Name for file with unmatched sequences [default:
            %(default)s]""", default="unmatched.sff")

    return parser


def main(parsed):
    """
    Command-line functionality

    * parse arguments
    * read barcodes
    * run
    """
    with parsed.barcode_file:
        barcodes = load_barcodes(parsed.barcode_file)

    with parsed.input_file:
        sequences = SeqIO.parse(parsed.input_file, 'sff')
        splitter = SequenceSplitter(parsed.primer, barcodes,
                                    parsed.make_dirs, parsed.unmatched_name)
        try:
            splitter.open()
            splitter.run(sequences)
        finally:
            splitter.close()
