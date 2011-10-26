"""
Helpers for running AmpliconNoise binaries
"""

import logging
import multiprocessing
import os
import shlex
import shutil
import subprocess
import tempfile

import argparse

def executable_transformer(use_m=False):
    """
    Append "M" to executable name if use_m = True
    """
    def m(executable):
        if use_m:
            return executable + 'M'
        else:
            return executable
    return m

def base_parser():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--stub', default=None,
            help="Stub for the output [default: base name of input file]")
    parser.add_argument('--temp-dir', default=tempfile.gettempdir(),
            help="Temporary directory [default: %(default)s]")
    parser.add_argument('--mpi-args', type=shlex.split,
            help="Arguments to pass to mpirun. [default: %(default)s]",
            default=['--np', multiprocessing.cpu_count()])
    parser.add_argument('-M', '--use-m', help="""Use M flavor of AmpliconNoise
            binaries [default: %(default)s]""", default=False,
            action='store_true')
    return parser


class NoiseRunner(object):
    """
    Runner for AmpliconNoise pipelines

    Runs commands in temporary directory, cleans on exit.
    """

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
        except Exception:
            logging.exception("Error fetching targets")
        finally:
            self._cleanup()
        return False
