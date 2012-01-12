import re
import sys

from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

from anoisetools.util import ambiguous_pattern


def build_parser(subparsers):
    """
    Adds command options to parser
    """
    parser = subparsers.add_parser('truncate',
             help="Remove tags and truncate length of FASTA files.",
             description="""Removes sequence <tag>,
trims remaining sequence to <length> from FASTA-formatted sequences passed
to stdin, printing to stdout.""")
    parser.add_argument('barcode', metavar='<tag>', help="""Sequence tag to
            remove if present. IUPAC ambiguous characters are accepted.""",
            type=ambiguous_pattern)
    parser.add_argument('length', metavar='<length>',
            help='Trim sequences to <length>', type=int)
    return parser

def trim(barcode, length, in_handle, out_handle):
    """
    Trim input sequences, write to out_handle

    Trims barcodes from FASTA-formatted sequences in in_handle, truncates
    sequences at provided length, writes to out_handle.
    """
    def inner(sequences, pattern):
        "Does the trimming"
        for sequence in sequences:
            s = str(sequence.seq)
            m = pattern.match(s)
            if m:
                yield sequence[m.end():length]
            else:
                print >> sys.stderr, "No match:", sequence.id
                yield sequence[:length]

    # Records are read TCAG, we trim the first 4 *flows*, which won't
    # detect duplicate G's.
    barcode = barcode.lstrip('G')

    pattern = re.compile('^{0}'.format(barcode))

    sequences = SeqIO.parse(in_handle, 'fasta')

    trimmed = inner(sequences, pattern)
    FastaWriter(out_handle, wrap=None).write_file(trimmed)

def main(parsed):
    """
    Calls ``trim(...)`` based on parsed arguments
    """
    trim(parsed.barcode, parsed.length, sys.stdin, sys.stdout)
