#!/usr/bin/env python

import logging
import os
import os.path
import multiprocessing
import shlex
import shutil
import subprocess
import tempfile

from anoisetools.scripts import clean, sff2raw
from anoisetools import anoiseio



def build_parser(subparsers):
    parser = subparsers.add_parser("pyronoise", help="Run PyroNoise")

    pnoise_opts = parser.add_argument_group("PyroNoise Options")
    pnoise_opts.add_argument('-s', '--sigma', type=float,
            default=60.0, help="""Initial cluster size (see Quince et al)
            [default: %(default)s]""")
    pnoise_opts.add_argument('-c', '--cutoff', type=float,
            default=0.01, help="""Initial cutoff (see Quince et al) [default:
            %(default)s]""")

    parser.add_argument('sff_file', help="""SFF file""")
    parser.add_argument('-v', '--verbose', default=logging.INFO,
            action="store_const", const=logging.DEBUG)

    parser.add_argument('--stub', default=None,
            help="Stub for the output [default: base name of SFF file]")
    parser.add_argument('--temp-dir', default=tempfile.gettempdir(),
            help="Temporary directory [default: %(default)s]")
    parser.add_argument('--mpi-args', type=shlex.split,
            help="Arguments to pass to mpirun. [default: %(default)s]",
            default=['--np', multiprocessing.cpu_count()])

    clean_group = parser.add_argument_group("Cleaning parameters")
    clean_group.add_argument('--min-flows', default=320, type=int,
            help="""Minimum flow count [default: %(default)s]""")
    clean_group.add_argument('--max-flows', default=720, type=int,
            help="""Maximum flow count [default: %(default)s]""")
    clean_group.add_argument('--max-empty', default=0, type=int,
            help="""Maximum empty flowgrams to accept before truncation
            [default: %(default)s]""")
    clean_group.add_argument('--primer', default='.',
            help="Regular expression to identify primer")

    return parser

def main(arguments):
    if not arguments.stub:
        arguments.stub = os.path.basename(
                os.path.splitext(arguments.sff_file)[0])

    pnoise_stub = arguments.stub + '-pnoise'
    targets = [pnoise_stub, pnoise_stub + '_cd.fa', pnoise_stub + '.mapping']


    runner = NoiseRunner(targets, temp_base=arguments.temp_dir,
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
                clean.invoke(reader, fa_fp, clean_dat, arguments.primer,
                        arguments.min_flows, arguments.max_flows,
                        arguments.max_empty)

        logging.info("Starting PyroNoise")
        run_pyronoise(runner, clean_dat, arguments.sigma,
                arguments.cutoff, arguments.stub)


    #snoise_results = run_seqnoise(pnoise_result['denoised_sequences'],
            #pnoise_result['mapping'], arguments.s_cluster_size, arguments.stub
            #+ '-snoise', arguments.np, arguments.nodes_file, tmp_dir)


class NoiseRunner(object):

    def __init__(self, target_files, target_dir=os.getcwd(), temp_base=None,
            mpi_cmd='mpirun', mpi_flags=None, cleanup=True):
        self.mpi_cmd = mpi_cmd
        self.mpi_flags = mpi_flags or []
        self.temp_base = temp_base
        self.target_files = target_files
        self.target_dir = target_dir
        self.temp_dir = None
        self.cleanup = cleanup

    def _mpi_command(self, cmd):
        return [self.mpi_cmd] + self.mpi_flags + cmd

    def path_join(self, fname):
        """
        Returns the full path to file identified by fname in the working temp
        dir.
        """
        return os.path.join(self.temp_dir, fname)

    def run(self, command, mpi=True, suppress_stdout=True,
            suppress_stderr=True, **kwargs):
        """
        Run the specified command
        """

        if not (self.temp_dir and os.path.isdir(self.temp_dir)):
            raise ValueError("Missing temporary directory: "
                    "{0}".format(self.temp_dir))

        if mpi:
            command = self._mpi_command(command)

        outputs = [i for i in ('stderr', 'stdout') if i not in kwargs]
        for i in ('stderr', 'stdout'):
            if not i in kwargs:
                kwargs[i] = tempfile.SpooledTemporaryFile()

        # Convert args to strings
        command = map(str, command)

        logging.info("Running: %s", " ".join(command))
        try:
            p = subprocess.Popen(command, cwd=self.temp_dir, **kwargs)
            p.wait()
            if p.returncode:
                raise subprocess.CalledProcessError("Error", p.returncode)
        except subprocess.CalledProcessError, e:
            if 'stderr' in outputs:
                kwargs['stderr'].seek(0)
                stderr = kwargs['stderr'].read()
            else:
                stderr = ''
            if 'stdout' in outputs:
                kwargs['stdout'].seek(0)
                stdout = kwargs['stdout'].read()
            else:
                stdout = ''

            logging.error("Error running '%s': returned with %s.\nstdout:\n%s\nstderr:%s",
                    ' '.join(command), p.returncode, stdout, stderr)
            raise SystemExit(1)

    def _setup(self):
        # Generate a temporary directory
        self.temp_dir = tempfile.mkdtemp(prefix='noise-', dir=self.temp_base)
        logging.info("Working in %s", self.temp_dir)

    def _fetch_targets(self):
        for f in self.target_files:
            logging.info("Moving %s to %s", f, self.target_dir)
            shutil.move(self.path_join(f), self.target_dir)

    def _cleanup(self):
        logging.info("Cleaning up")
        # Move all the results in
        if self.temp_dir:
            if self.cleanup:
                logging.info("Removing %s", self.temp_dir)
                shutil.rmtree(self.temp_dir)
            else:
                logging.info("Keeping %s", self.temp_dir)
        self.temp_dir = None

    def __enter__(self):
        self._setup()
        return self

    def __exit__(self, err_type, err_value, err_traceback):
        try:
            if not err_type:
                self._fetch_targets()
        except Exception, e:
            logging.exception("Error fetching targets")
        finally:
            self._cleanup()
        return False


def run_pyronoise(runner, dat_path, sigma, cutoff, output_stub):
    pnoise_stub = output_stub + '-pnoise'
    runner.run(['PyroDist', '-in', dat_path, '-out', output_stub])
    runner.run(['FCluster', '-in', output_stub + '.fdist',
        '-out', output_stub + '-initial'])
    runner.run(['PyroNoise', '-din', dat_path, '-out', pnoise_stub,
            '-lin', output_stub + '-initial.list', '-s', sigma,
            '-c', cutoff])

