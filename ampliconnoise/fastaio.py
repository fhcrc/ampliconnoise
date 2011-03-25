import collections
import itertools
import textwrap

Sequence = collections.namedtuple('Sequence', ('id', 'seq'))


def parse_fasta(fp):
    """
    Parses a fasta file, generating Sequences
    """
    def is_header(line):
        return line and line.startswith(">")

    # Generate
    lines = (line.strip().rstrip('\n') for line in fp)

    # Grab the first line
    while True:
        line = next(lines)
        if line.startswith('>'):
            break
        else:
            raise ValueError("Unexpected line: {0}".format(line))

    while True:
        header = line
        seq_lines = []
        while True:
            line = next(lines)
            if line.startswith('>'):
                break
            else:
                seq_lines.append(line)
        yield Sequence(header, ''.join(seq_lines))


def write_fasta(sequences, fp, wrap=None):
    """
    Write sequences to fp.

    If wrap is given, sequences will be wrapped to lines of length wrap.
    """
    count = 0
    for sequence in sequences:
        header = '>' + sequence.id
        seq = sequence.seq
        if wrap is not None:
            if wrap is True:
                wrap = 80
            seq = '\n'.join(textwrap.wrap(sequence, wrap))
        print >> fp, header
        fp.write(seq)
        count += 1

    return count

