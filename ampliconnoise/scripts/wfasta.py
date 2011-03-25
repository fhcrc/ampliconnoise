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
_HEADER_RE = re.compile(r'^>(.*?)_(\d+)_(\d+)$')


def build_parser(parent_subparsers):
    parser = parent_subparsers.add_parser('wfasta',
            help='Work with weighted FASTA output by SeqNoise')
    subparsers = parser.add_subparsers()
    rep_parser = subparsers.add_parser('repeat',
            help='Repeat sequences <weight> times')
    rep_parser.add_argument('--wrap', default=None, type=int,
            help="Wrap sequences to <wrap> lines. Default: no wrap")
    rep_parser.set_defaults(subfunc=invoke_repeat)

    tab_parser = subparsers.add_parser('tabular',
            help='Convert to tab-delimited')
    tab_parser.set_defaults(subfunc=invoke_tabular)

    parser.add_argument('--min-frequency', type=int,
        help="Minimum frequency for output")
    parser.add_argument('infile', type=argparse.FileType('r'),
            help="Infile (default: stdin)", default=sys.stdin)
    parser.add_argument('outfile', type=argparse.FileType('w'),
            help="Outfile (default: stdout)", default=sys.stdout)

    return parser

def main(parsed):
    parsed.subfunc(parsed)

def _parse_wfasta_header(line):
    m = _HEADER_RE.match(line)
    if not m:
        raise ValueError("Invalid wfasta header: {0}".format(line))
    i, index, freq = m.groups()
    return WeightedFastaHeader(i, int(index), int(freq))


def repeat(sequences, min_frequency=None):
    for record in sequences:
        parsed_header = _parse_wfasta_header(record.id)
        # Skip if the expected frequencies are not present
        if min_frequency and parsed_header.frequency < min_frequency:
            continue
        for i in xrange(parsed_header.frequency):
            new_header = '{0}_{1}'.format(record.id, i)
            yield fastaio.Sequence(new_header, record.seq)


def invoke_repeat(parsed_args):
    """
    Invokes repeat command
    """
    with parsed_args.infile:
        with parsed_args.outfile:
            reader = fastaio.parse_fasta(parsed_args.infile)
            repeated_sequences = repeat(reader, parsed_args.min_frequency)
            fastaio.write_fasta(repeated_sequences, parsed_args.outfile,
                                parsed_args.wrap)

def invoke_tabular(parsed_args):
    with parsed_args.infile:
        with parsed_args.outfile:
            writer = csv.writer(parsed_args.outfile, delimiter='\t',
                    lineterminator='\n', quoting=csv.QUOTE_NONE)
            writer.writerow(('id', 'index', 'frequency', 'sequence'))
            reader = fastaio.parse_fasta(parsed_args.infile)
            for sequence in reader:
                parsed_header = _parse_wfasta_header(sequence.id)
                if (parsed_args.min_frequency is None or
                        parsed_header.frequency >= parsed_args.min_frequency):
                    writer.writerow((parsed_header.id, parsed_header.index,
                                     parsed_header.frequency, sequence.seq))

