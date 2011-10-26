#!/usr/bin/env python

import logging
import os
import os.path

from anoisetools.scripts import clean, sff2raw
from anoisetools import anoiseio, run, util


def build_parser(subparsers):
    parser = subparsers.add_parser("pyronoise", help="Run PyroNoise",
            parents=[run.base_parser()])

    pnoise_opts = parser.add_argument_group("PyroNoise Options")
    pnoise_opts.add_argument('-s', '--sigma', type=float,
            default=60.0, help="""Initial cluster size (see Quince et al)
            [default: %(default)s]""")
    pnoise_opts.add_argument('-c', '--cutoff', type=float,
            default=0.01, help="""Initial cutoff (see Quince et al) [default:
            %(default)s]""")
    pnoise_opts.add_argument('-rin', '--lookup-file',
            help="""Lookup file""")

    parser.add_argument('sff_file', help="""SFF file""")
    parser.add_argument('-v', '--verbose', default=logging.INFO,
            action="store_const", const=logging.DEBUG)

    clean_group = parser.add_argument_group("Cleaning parameters")
    clean_group.add_argument('--min-flows', default=320, type=int,
            help="""Minimum flow count [default: %(default)s]""")
    clean_group.add_argument('--max-flows', default=720, type=int,
            help="""Maximum flow count [default: %(default)s]""")
    clean_group.add_argument('--max-empty', default=0, type=int,
            help="""Maximum empty flowgrams to accept before truncation
            [default: %(default)s]""")
    clean_group.add_argument('--primer', default='.',
            help="Primer used", type=util.ambiguous_pattern)

    return parser

def main(arguments):
    if not arguments.stub:
        arguments.stub = os.path.basename(
                os.path.splitext(arguments.sff_file)[0])

    pnoise_stub = arguments.stub + '-pnoise'
    targets = [pnoise_stub, pnoise_stub + '_cd.fa', pnoise_stub + '.mapping']

    runner = run.NoiseRunner(targets, temp_base=arguments.temp_dir,
        mpi_flags=arguments.mpi_args)

    with runner:
        # Clean
        logging.info("Extracting clean flows")
        raw_file = runner.path_join(arguments.stub + '.raw')
        with open(raw_file, 'w') as fp:
            sff2raw.sff_to_raw(arguments.sff_file, fp)
        clean_dat = runner.path_join(arguments.stub + '.dat')
        with open(raw_file) as raw_fp:
            reader = anoiseio.AnoiseRawReader(raw_fp)
            with open(runner.path_join(arguments.stub + '.fa'), 'w') as fa_fp:
                passed, failed = clean.invoke(reader, fa_fp, clean_dat,
                        arguments.primer, arguments.min_flows,
                        arguments.max_flows, arguments.max_empty)

                if not passed:
                    raise AssertionError(
                            "No sequences passed flowgram cleaning!")

        logging.info("Starting PyroNoise")
        run_pyronoise(runner, clean_dat, arguments.sigma,
                arguments.cutoff, arguments.stub, arguments.use_m,
                lookup_name=arguments.lookup_file)


def run_pyronoise(runner, dat_path, sigma, cutoff, output_stub, use_m=False,
        lookup_name=None):
    m = run.executable_transformer(use_m)
    pnoise_stub = output_stub + '-pnoise'
    runner.run(['PyroDist', '-in', dat_path, '-out', output_stub])
    runner.run([m('FCluster'), '-in', output_stub + '.fdist',
        '-out', output_stub + '-initial'])
    cmd = [m('PyroNoise'),
           '-din', dat_path,
           '-out', pnoise_stub,
           '-lin', output_stub + '-initial.list',
           '-s', sigma,
           '-c', cutoff]
    if lookup_name:
        cmd.extend(('-rin', lookup_name))
    runner.run(cmd)
