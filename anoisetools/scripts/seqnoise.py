import logging
import multiprocessing
import tempfile
import shlex
import os
import os.path

# TODO: Refactor
from anoisetools.scripts.pyronoise import NoiseRunner


def build_parser(subparsers):
    logging.basicConfig(level=logging.INFO)
    parser = subparsers.add_parser("seqnoise", help="Run SeqNoise")

    pnoise_opts = parser.add_argument_group("SeqNoise Options")
    pnoise_opts.add_argument('-s', '--sigma', type=float,
            default=30.0, help="""Sigma (see Quince et al)
            [default: %(default)s]""")
    pnoise_opts.add_argument('-c', '--cutoff', type=float,
            default=0.08, help="""Initial cutoff (see Quince et al) [default:
            %(default)s]""")

    parser.add_argument('-v', '--verbose', default=logging.INFO,
        action="store_const", const=logging.DEBUG)

    parser.add_argument('fasta_file', help="""PyroNoised fasta file""")
    parser.add_argument('mapping_file', help="""PyroNoise .mapping file""")

    parser.add_argument('--stub', default=None,
            help="Stub for the output [default: base name of Fasta file]")
    parser.add_argument('--temp-dir', default=tempfile.gettempdir(),
            help="Temporary directory [default: %(default)s]")
    parser.add_argument('--mpi-args', type=shlex.split,
            help="Arguments to pass to mpirun. [default: %(default)s]",
            default=['--np', multiprocessing.cpu_count()])

    return parser

def main(arguments):
    if not arguments.stub:
        arguments.stub = os.path.basename(
                os.path.splitext(arguments.fasta_file)[0])
    if arguments.stub.endswith('-pnoise'):
        arguments.stub = arguments.stub[:-7]

    snoise_stub = arguments.stub + '-snoise'
    targets = [snoise_stub, snoise_stub + '_cd.fa', snoise_stub + '.mapping']

    runner = NoiseRunner(targets, temp_base=arguments.temp_dir,
        mpi_flags=arguments.mpi_args, write_stdout=arguments.verbose)

    with runner:
        logging.info("Running SeqNoise")
        run_seqnoise(runner, arguments.fasta_file, arguments.mapping_file,
            arguments.sigma, arguments.cutoff, arguments.stub + '-snoise')

def run_seqnoise(runner, fasta_file, pnoise_mapping, sigma, cutoff, stub):
    # SeqDist
    seqdist_out = runner.path_join(stub + '.seqdist')
    with open(seqdist_out, 'w') as fp:
        runner.run(['SeqDist', '-in', os.path.abspath(fasta_file)], stdout=fp)

    runner.run(['FCluster', '-in', seqdist_out, '-out', stub])
    fcluster_list = stub + '.list'

    # SeqNoise
    runner.run(['SeqNoise',
                '-in', os.path.abspath(fasta_file),
                '-din', seqdist_out,
                '-out', stub,
                '-lin', fcluster_list,
                '-min', os.path.abspath(pnoise_mapping),
                '-c', cutoff,
                '-s', sigma])

