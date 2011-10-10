"""
Split input FASTA files
"""

import argparse
import collections
import csv
import errno
import logging
import re
import multiprocessing
import os
import os.path
import Queue

from Bio import SeqIO


# Map from Ambiguous Base to regex
_AMBIGUOUS_MAP = {
       'R': '[GA]',
       'Y': '[TC]',
       'K': '[GT]',
       'M': '[AC]',
       'S': '[GC]',
       'W': '[AT]',
       'B': '[GTC]',
       'D': '[GAT]',
       'H': '[ACT]',
       'V': '[GCA]',
       'N': '[AGCT]',
}

class SequenceWriter(object):
    """
    SequenceWriter writes sequences.

    Run as a subprocess, instances are initialized with a
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

                # None sent as a sentinal value for completed job
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

def _ambiguous_regex(sequence_str):
    return re.compile(''.join(_AMBIGUOUS_MAP.get(c, c) for c in sequence_str))

def load_barcodes(fp):
    d = collections.defaultdict(list)
    reader = csv.reader(fp)
    for record in reader:
        name, barcode, primer = record[:3]
        if barcode in d:
            raise ValueError("Duplicate key: {0}".format(barcode))
        if name in (name for p, name in d.values()):
            raise ValueError("Duplicate value: {0}".format(name))
        d[barcode].append((_ambiguous_regex(primer), name))

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

    def __init__(self, barcodes, make_dirs=False,
            unmatched_name="unknown.sff", output_format='sff'):
        self.barcodes = barcodes

        min_len = min(len(barcode) for barcode in barcodes)
        max_len = max(len(barcode) for barcode in barcodes)
        if min_len != max_len:
            raise ValueError("Unequal barcode lengths are not supported")

        self.barcode_length = max_len
        self.output_format = output_format
        self.make_dirs = make_dirs
        self.unmatched_name = unmatched_name

        self._output_files = collections.defaultdict(dict)
        self._barcode_primer_name = {}

    def open(self):
        """
        Create a subprocess running a SequenceWriter for each barcode.
        """
        if self._output_files:
            raise AlreadyOpenError("SequenceSplitter has already been opened.")

        writer_templates = []

        # Build a list of barcode, filepath tuples to construct writers from
        for barcode, primer_list in self.barcodes.items():
            for primer, identifier in primer_list:
                path = identifier + '.' + self.output_format
                if self.make_dirs:
                    path = os.path.join(identifier, path)
                try:
                    os.makedirs(os.path.dirname(path))
                except OSError, e:
                    if e.errno != errno.EEXIST:
                        raise

                self._barcode_primer_name[barcode, primer] = identifier

            writer_templates.append((barcode, primer, path))

        # Special case: unmatched barcodes
        writer_templates.append((None, self.unmatched_name))

        for barcode, primer, path in writer_templates:
            queue = multiprocessing.Queue()
            writer = SequenceWriter(queue, path, self.output_format)

            # Create a process to write output records
            process = multiprocessing.Process(target=writer.write)
            process.start()
            self._output_files[barcode][primer] = queue, writer, process

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        """
        Close all the output files
        """
        for queue, writer, process in (v for i in self._output_files.values()
                                       for v in i.values()):
            # Sentinal - no more records to process
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
        # TODO: Count based on name
        counter = collections.Counter()
        for sequence in sequences:
            seq_str = str(sequence.seq)[self.barcode_start:]
            barcode = seq_str[:self.barcode_length]

            queue = self._output_files[None][0]
            if barcode in self.barcodes:
                primers = self.primers[barcode]
                for primer_re in primers:
                    if primer_re.match(seq_str[self.barcode_length:]):
                        counter[barcode] += 1
                        queue = self._output_files.get(
                                barcode, self._output_files[None])[0]


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
    parser = subparsers.add_parser('split', help="Split an .sff")

    parser.add_argument("barcode_file", type=argparse.FileType('r'),
            help="""Path to barcode file with barcode_id,barcode_seq,primer
            tuples. Ambiguous characters are accepted.""")
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
        splitter = SequenceSplitter(barcodes,
                                    parsed.make_dirs, parsed.unmatched_name)
        with splitter:
            splitter.run(sequences)
