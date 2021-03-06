import logging
import os
import os.path

from anoisetools import run


def build_parser(subparsers):
    parser = subparsers.add_parser("seqnoise", help="Run SeqNoise",
            parents=[run.base_parser()])

    snoise_opts = parser.add_argument_group("SeqNoise Options")
    snoise_opts.add_argument('-s', '--sigma', type=float,
            default=30.0, help="""Sigma (see Quince et al)
            [default: %(default)s]""")
    snoise_opts.add_argument('-c', '--cutoff', type=float,
            default=0.08, help="""Initial cutoff (see Quince et al) [default:
            %(default)s]""")
    snoise_opts.add_argument('-rin', '--lookup-file',
            help="""Lookup file""")

    parser.add_argument('-v', '--verbose', default=logging.INFO,
        action="store_const", const=logging.DEBUG)

    parser.add_argument('fasta_file', help="""PyroNoised fasta file""")
    parser.add_argument('mapping_file', help="""PyroNoise .mapping file""")

    return parser

def main(arguments):
    if not arguments.stub:
        arguments.stub = os.path.basename(
                os.path.splitext(arguments.fasta_file)[0])
    if arguments.stub.endswith('-pnoise_cd'):
        arguments.stub = arguments.stub[:-10]

    snoise_stub = arguments.stub + '-snoise'
    targets = [snoise_stub + '_cd.fa', snoise_stub + '.mapping']

    runner = run.NoiseRunner(targets, temp_base=arguments.temp_dir,
        mpi_flags=arguments.mpi_args)

    with runner:
        logging.info("Running SeqNoise")
        run_seqnoise(runner, arguments.fasta_file, arguments.mapping_file,
            arguments.sigma, arguments.cutoff, arguments.stub + '-snoise',
            use_m=arguments.use_m, lookup_file=arguments.lookup_file)

def run_seqnoise(runner, fasta_file, pnoise_mapping, sigma, cutoff, stub,
        use_m=False, lookup_file=None):
    m = run.executable_transformer(use_m)
    # SeqDist
    seqdist_out = runner.path_join(stub + '.seqdist')
    with open(seqdist_out, 'w') as fp:
        runner.run(['SeqDist', '-in', os.path.abspath(fasta_file)], stdout=fp)

    runner.run([m('FCluster'), '-in', seqdist_out, '-out', stub])
    fcluster_list = stub + '.list'

    # SeqNoise
    cmd = [m('SeqNoise'),
                '-in', os.path.abspath(fasta_file),
                '-din', seqdist_out,
                '-out', stub,
                '-lin', fcluster_list,
                '-min', os.path.abspath(pnoise_mapping),
                '-c', cutoff,
                '-s', sigma]
    if lookup_file:
        cmd.extend(('-rin', lookup_file))
    runner.run(cmd)
