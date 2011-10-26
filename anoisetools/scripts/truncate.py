import re
import sys

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

    Note:
    Sequences passed must be on a single line.
    """
    # Records are read TCAG, we trim the first 4 *flows*, which won't
    # detect duplicate G's.
    barcode = barcode.lstrip('G')
    lines = (line.rstrip() for line in in_handle)

    pattern = re.compile('{0}(.*)'.format(barcode))



    for line in lines:
        if line.startswith('>'):
            # Fasta header
            print >> out_handle, line
        else:
            m = pattern.match(line)

            if m:
                line = m.group(1)

            print >> out_handle, line[:length]


def main(parsed):
    """
    Calls ``trim(...)`` based on parsed arguments
    """
    trim(parsed.barcode, parsed.length, sys.stdin, sys.stdout)
