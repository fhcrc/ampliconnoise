"""
Tools for working with weighted FASTA output from SeqNoise
"""
import csv
import argparse
import collections
import re
import sys

from ampliconnoise import fastaio


WeightedFastaHeader = collections.namedtuple('WeightedFastaHeader',
        ('id', 'index', 'frequency'))
_HEADER_RE = re.compile(r'^(.*?)_(\d+)_(\d+)$')


class _WFastaWriter(object):

    def __init__(self, fp, min_frequency, repeat):
        self.fp = fp
        self.min_frequency = min_frequency
        self.repeat = repeat

    def writerecord(self, record):
        header = _parse_wfasta_header(record.id)
        if (self.min_frequency is None or
                header.frequency >= self.min_frequency):
            self.write(record)


class _WFastaTabularWriter(_WFastaWriter):

    def __init__(self, fp, min_frequency, repeat):
        super(_WFastaTabularWriter, self).__init__(fp, min_frequency, repeat)
        self.writer = csv.writer(self.fp, delimiter='\t', lineterminator='\n')
        self.writer.writerow(('id', 'index', 'frequency', 'sequence'))

    def write(self, sequence):
        parsed_header = _parse_wfasta_header(sequence.id)
        self.writer.writerow((parsed_header.id, parsed_header.index,
                              parsed_header.frequency, sequence.seq))


class _WFastaFastaWriter(_WFastaWriter):

    def write(self, sequence):
        parsed_header = _parse_wfasta_header(sequence.id)
        seqid = sequence.id

        towrite = []
        if self.repeat is True:
            for i in xrange(parsed_header.frequency):
                new_header = '{0}_{1}'.format(record.id, i)
                towrite.append(fastaio.Sequence(new_header, record.seq))
        else:
            towrite.append(sequence)

        fastaio.write_fasta(towrite, self.fp)


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
    parser.add_argument('infile', type=argparse.FileType('r'),
            help="Infile (default: stdin)", default=sys.stdin)
    parser.add_argument('outfile', type=argparse.FileType('w'),
            help="Outfile (default: stdout)", default=sys.stdout)

    return parser


def main(parsed):
    with parsed.infile:
        with parsed.outfile:
            writer_class = _OUTPUT_FORMATS[parsed.output_format]
            writer = writer_class(parsed.outfile, parsed.min_frequency,
                                  parsed.repeat)

            for sequence in fastaio.parse_fasta(parsed.infile):
                writer.writerecord(sequence)


def _parse_wfasta_header(line):
    """
    Parse wfasta header from the fasta id
    """
    m = _HEADER_RE.match(line)
    if not m:
        raise ValueError("Invalid wfasta header: {0}".format(line))
    i, index, freq = m.groups()
    return WeightedFastaHeader(i, int(index), int(freq))
