"""
Tools for working with weighted FASTA output from SeqNoise
"""
import csv
import argparse
import collections
import re
import sys

from Bio import SeqIO


WeightedFastaHeader = collections.namedtuple('WeightedFastaHeader',
        ('id', 'index', 'frequency'))

_HEADER_RE = re.compile(r'^(.*?)_(\d+)_(\d+)$')

def _parse_wfasta_header(line):
    """
    Parse wfasta header from the fasta id field
    """
    m = _HEADER_RE.match(line)
    if not m:
        raise ValueError("Invalid wfasta header: {0}".format(line))
    i, index, freq = m.groups()
    return WeightedFastaHeader(i, int(index), int(freq))

class _WFastaWriter(object):
    """
    base class for WFasta writers
    """

    def __init__(self, fp, min_frequency=None, repeat=False):
        """
        Initialize the writer with a destination file pointer,
        minimum frequency to accept, and a boolean indicating whether
        to repeat the sequence.
        """
        self.fp = fp
        self.min_frequency = min_frequency
        self.repeat = repeat

    def writerecord(self, record):
        """
        Write a record to the output file.

        This method should not need to be overridden - the real action
        is in _WFastaWriter.write.
        """
        header = _parse_wfasta_header(record.id)
        if (self.min_frequency is None or
                header.frequency >= self.min_frequency):
            self.write(record)

    def write(self, sequence):
        """
        Writes a sequence object to the output file.
        """
        raise NotImplementedError("Override in subclass")


class _WFastaTabularWriter(_WFastaWriter):
    """
    Writer for tab-delimited text, with a header
    """

    def __init__(self, fp, min_frequency, repeat):
        super(_WFastaTabularWriter, self).__init__(fp, min_frequency, repeat)
        self.writer = csv.writer(self.fp, delimiter='\t', lineterminator='\n')
        self.writer.writerow(('id', 'index', 'frequency', 'sequence'))

    def write(self, sequence):
        """
        Writes a sequence object to the output file in tab-delimited format.
        """
        parsed_header = _parse_wfasta_header(sequence.id)
        self.writer.writerow((parsed_header.id, parsed_header.index,
                              parsed_header.frequency, sequence.seq))


class _WFastaFastaWriter(_WFastaWriter):
    """
    Write the output in FASTA.
    """

    def write(self, sequence):
        parsed_header = _parse_wfasta_header(sequence.id)
        seqid = sequence.id

        if self.repeat is True:
            current_id = record.id
            for i in xrange(parsed_header.frequency):
                new_header = '{0}_{1}'.format(current_id, i)
                sequence.id = new_header
                SeqIO.write([sequence], fp, 'fasta')
        else:
            SeqIO.write([sequence], fp, 'fasta')


_OUTPUT_FORMATS = {'fasta': _WFastaFastaWriter,
                   'tabular': _WFastaTabularWriter}


def build_parser(parent_subparsers):
    """
    Build and return a command parser
    """
    parser = parent_subparsers.add_parser('wfasta',
            help='Work with weighted FASTA output by SeqNoise')
    parser.add_argument('--min-frequency', type=int,
        help="Minimum frequency for output")
    parser.add_argument('--output-format', choices=_OUTPUT_FORMATS.keys(),
            help="Output format (choices: %(choices)s, default: %(default)s)",
            default="fasta")
    parser.add_argument('--repeat', action="store_true", default=False,
            help="Repeat each record [frequency] times (default %(default)s)")
    parser.add_argument('--infile', type=argparse.FileType('r'),
            help="Infile (default: stdin)", default=sys.stdin)
    parser.add_argument('--outfile', type=argparse.FileType('w'),
            help="Outfile (default: stdout)", default=sys.stdout)

    return parser


def main(parsed):
    """
    Filters / converts the input file based on parsed command line arguments
    """
    with parsed.infile:
        with parsed.outfile:
            writer_class = _OUTPUT_FORMATS[parsed.output_format]
            writer = writer_class(parsed.outfile, parsed.min_frequency,
                                  parsed.repeat)

            for sequence in SeqIO.parse(parsed.infile, 'fasta'):
                writer.writerecord(sequence)
