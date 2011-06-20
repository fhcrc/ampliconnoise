import argparse
import glob
import os.path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def build_parser(subparsers):
    parser = subparsers.add_parser("redup",
            help="Reduplicate post-SeqNoise sequences")
    parser.add_argument('fasta',
            type=argparse.FileType('r'), help="""Post-SeqNoise FASTA file""")
    parser.add_argument('pnoise_mapping_dir',
            help="PyroNoise output dir (with i_*fa files)")
    parser.add_argument('snoise_mapping_dir',
            help="SeqNoise -> PyroNoise ID mapping")
    parser.add_argument('output', type=argparse.FileType('w'))

    return parser


def parse_seqnoise_mapping(fp):
    """
    Parse an AmpliconNoise mapping file, consisting of:

        ID<comma-separated original IDs>

    returns a dictionary mapping from ID to a list of original IDs
    """
    lines = (line.rstrip('\n').split(None, 1) for line in fp)
    result = ((new_id, collapsed_ids.split(','))
              for new_id, collapsed_ids in lines)
    return dict(result)


def parse_noise_mapping(pnoise_dir):
    """
    Parse an AmpliconNoise mapping file, consisting of:

        ID <comma-separated original IDs>

    returns a dictionary
    """
    result = {}
    for f in glob.glob(os.path.join(pnoise_dir, 'i_[0-9]*.fa')):
        with open(f) as fp:
            ids = [s.id for s in SeqIO.parse(fp, 'fasta')]
        result[ids[0].lstrip('>')] = ids[1:]

    return result


def redup(input_sequences, mapping):
    for s in input_sequences:
        mapped = mapping[s.id]
        for m in mapped:
            yield SeqRecord(s.seq, id=m, description="")


def main(parsed):
    # Get mappings
    snoise_mapping = parse_noise_mapping(parsed.snoise_mapping_dir)
    pnoise_mapping = parse_noise_mapping(parsed.pnoise_mapping_dir)

    # Merge, mapping from snoise to original sequence
    mapping = {}
    for snoise_id, pnoise_ids in snoise_mapping.items():
        print snoise_id
        for p in pnoise_ids:
            # Extract ID from sequences matching p5z1r04-pnoise_30_7
            try:
                mapping[snoise_id].extend(pnoise_mapping[p])
            except KeyError:
                mapping[snoise_id] = pnoise_mapping[p]

    with parsed.fasta:
        sequences = SeqIO.parse(parsed.fasta, 'fasta')
        with parsed.output:
            SeqIO.write(redup(sequences, mapping), parsed.output, 'fasta')
