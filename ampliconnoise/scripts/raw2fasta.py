"""
Convert from .raw to .fasta files
"""
import argparse
import os.path
import sys

from ampliconnoise import anoiseio


def build_parser(subparsers):
    """
    Add arguments
    """
    parser = subparsers.add_parser('raw2fasta', help="""Convert from .raw to
.fasta""")
    parser.add_argument('infile', metavar='<raw file>',
            type=argparse.FileType('r'), help="Input file")
    return parser


def raw_to_fasta(infile, outfile):
    """
    Reads raw records from infile, converts to fasta, writes to outfile
    """
    reader = anoiseio.AnoiseRawReader(infile)
    for record in reader:
        print >> outfile, '>{0}'.format(record.identifier)
        print >> outfile, record.bases_from_flows()


def main(parsed):
    """
    Runs raw_to_fasta with given arguments
    """
    noext = os.path.splitext(parsed.infile.name)[0]
    outfile_path = noext + '.fasta'
    with open(outfile_path, 'w') as outfile:
        with parsed.infile:
            raw_to_fasta(parsed.infile, outfile)
