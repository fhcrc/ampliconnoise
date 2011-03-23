import argparse
import re
import sys


def trim(barcode, length, in_handle, out_handle):
    # trim leading G's #TODO: why
    # barcode = barcode.lstrip('G')
    lines = (line.rstrip() for line in in_handle)
    pattern = re.compile('{0}(.*)'.format(barcode))

    for line in lines:
        if line.startswith('>'):
            # Fasta header
            print >>out_handle, line
        else:
            m = pattern.match(line)

            if not m:
                raise ValueError("Unknown sequence: {0}".format(line))

            seq = m.group(1)
            print >> out_handle, seq[:length]

def main(args=sys.argv[1:]):

    parser = argparse.ArgumentParser(description="""Removes sequence <tag>,
trims remaining sequence to <length> from FASTA-formatted sequences passed
to stdin, printing to stdout.""")
    parser.add_argument('barcode', metavar='<tag>',
            help='Sequence tag')
    parser.add_argument('length', metavar='<length>',
            help='Trim sequences to <length>', type=int)

    parsed = parser.parse_args(args)

    trim(parsed.barcode, parsed.length, sys.stdin, sys.stdout)
