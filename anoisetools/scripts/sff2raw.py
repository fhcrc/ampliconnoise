"""
Convert from .raw to .fasta files
"""
import argparse
import os.path

from anoisetools import anoiseio

from Bio import SeqIO


def build_parser(subparsers):
    """
    Add arguments
    """
    parser = subparsers.add_parser('sff2raw', help="""Convert from .sff to
    .raw""")
    parser.add_argument('infile', metavar='SFF',
            type=argparse.FileType('rb'), help="Input file")
    parser.add_argument('outfile', metavar='rawfile',
            type=argparse.FileType('w'), help="output .raw", nargs='?')
    return parser


def sff_to_raw(infile, outfile, identifier=None):
    """
    Reads raw records from infile, converts to fasta, writes to outfile
    """
    sequences = SeqIO.parse(infile, 'sff')
    writer = anoiseio.AnoiseRawWriter(outfile, identifier or '454Reads')
    try:
        writer.write_records(sequences)
    finally:
        writer.close()


def main(parsed):
    """
    Runs sff_to_raw with given arguments
    """
    if not parsed.outfile:
        noext = os.path.splitext(parsed.infile.name)[0]
        outfile_path = noext + '.raw'
        parsed.outfile = open(outfile_path, 'w')
    with parsed.infile:
        with parsed.outfile:
            sff_to_raw(parsed.infile, parsed.outfile)
