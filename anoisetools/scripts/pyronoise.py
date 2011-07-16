#!/usr/bin/env python

import logging
import os
import os.path
import multiprocessing
import shlex
import shutil
import subprocess
import tempfile


def build_parser(subparsers):
    logging.basicConfig(level=logging.INFO)
    parser = subparsers.add_parser("pyronoise", help="Run PyroNoise")

    pnoise_opts = parser.add_argument_group("PyroNoise Options")
    pnoise_opts.add_argument('-s', '--p-cluster-size', type=float,
            default=60.0, help="""Initial cluster size (see Quince et al)
            [default: %(default)s]""")
    pnoise_opts.add_argument('-c', '--p-initial-cutoff', type=float,
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

    return parser

def main(arguments):


    # Override check_call
    sp_check_call = subprocess.check_call
    def check_call(cmd, *args, **kwargs):
        logging.info(' '.join(cmd))
        sp_check_call(cmd, *args, **kwargs)

    subprocess.check_call = check_call

    tmp_dir = create_temp_dir(arguments.temp_dir)
    logging.info("Working in temporary directory: %s", tmp_dir)

    if not arguments.stub:
        arguments.stub = os.path.basename(os.path.splitext(arguments.stub)[0])

    clean_dat = os.path.join(tmp_dir, arguments.stub + '.dat')
    logging.info("Writing clean flows to %s", clean_dat)
    clean_sff(arguments.sff_file, arguments.clean_params, clean_dat)

    pnoise_result = run_pyronoise(clean_dat, arguments.p_cluster_size,
            arguments.p_initial_cutoff, arguments.np, arguments.nodes_file,
            arguments.stub, cwd=tmp_dir)
    for v in pnoise_result.values():
        logging.info(v)
        shutil.move(v, os.getcwd())

    pnoise_result = {k:os.path.abspath(os.path.basename(v))
                     for k, v in pnoise_result.items()}

    # Clear
    logging.info("Removing %s", tmp_dir)
    shutil.rmtree(tmp_dir)

    tmp_dir = create_temp_dir()

    snoise_results = run_seqnoise(pnoise_result['denoised_sequences'],
            pnoise_result['mapping'], arguments.s_cluster_size, arguments.stub
            + '-snoise', arguments.np, arguments.nodes_file, tmp_dir)

    for v in snoise_results.values():
        logging.info(v)
        shutil.move(v, os.getcwd())

    logging.info("removing %s", tmp_dir)
    shutil.rmtree(tmp_dir)

def create_temp_dir(temp_dir=None):
    return tempfile.mkdtemp(prefix='anoise', dir=temp_dir)

def mpirun(command, np=1, nodes_file=None):
    cmd = ['mpirun', '-np', str(np)]
    if nodes_file:
        cmd.extend(['-machinefile', nodes_file])
    cmd.append(command)
    return cmd

def clean_sff(sff_path, clean_params, outfile):
    bn = os.path.splitext(outfile)[0]
    rawpath = bn + '.raw'
    cmd = ['sff2raw', os.path.basename(sff_path), sff_path, rawpath]
    subprocess.check_call(cmd)
    cmd = ['anoise', 'clean', '--input', rawpath, '.', bn]
    cmd.extend(clean_params)
    subprocess.check_call(cmd)

def run_pyronoise(dat_path, p_cluster_size, p_initial_cutoff,
        np, nodes_file, output_stub, cwd=None):

    # PyroDist
    pyrodist_cmd = mpirun('PyroDist', np, nodes_file)
    pyrodist_cmd.extend(['-in', dat_path, '-out', output_stub])
    subprocess.check_call(pyrodist_cmd, cwd=cwd)

    # FCluster
    fcluster_cmd = mpirun('FCluster', np, nodes_file)
    fcluster_cmd.extend(['-in', output_stub + '.fdist',
                         '-out', output_stub + '-initial'])
    subprocess.check_call(fcluster_cmd, cwd=cwd)

    # PyroNoise
    pnoise_stub = output_stub + '-pnoise'
    pyronoise_cmd = mpirun('PyroNoise', np, nodes_file)
    pyronoise_cmd.extend(['-din', dat_path, '-out', pnoise_stub,
            '-lin', output_stub + '-initial.list', '-s', p_cluster_size,
            '-c', p_initial_cutoff])
    subprocess.check_call(map(str, pyronoise_cmd), cwd=cwd)

    # Return output needed by results
    return {'sequence_map': os.path.join(cwd, pnoise_stub),
            'denoised_sequences': os.path.join(cwd, pnoise_stub + '_cd.fa'),
            'mapping': os.path.join(cwd, pnoise_stub + '.mapping')}

def run_seqnoise(fasta_file, pnoise_mapping, s_cluster_size, stub,
        np, nodes_file, cwd=None):

    # SeqDist
    seqdist_cmd = mpirun('SeqDist', np, nodes_file)
    seqdist_out = os.path.join(cwd or '.', stub + '.seqdist')
    seqdist_cmd.extend(['-in', fasta_file])
    with open(seqdist_out, 'w') as fp:
        subprocess.check_call(seqdist_cmd, stdout=fp, cwd=cwd)

    # FCluster
    fcluster_cmd = mpirun('FCluster', np, nodes_file)
    fcluster_cmd.extend(['-in', seqdist_out, '-out', stub])
    subprocess.check_call(fcluster_cmd, cwd=cwd)
    fcluster_list = stub + '.list'

    # SeqNoise
    seqnoise_cmd = mpirun('SeqNoise', np, nodes_file)
    seqnoise_cmd.extend(['-in', fasta_file,
                         '-din', seqdist_out,
                         '-out', stub,
                         '-lin', fcluster_list,
                         '-min', pnoise_mapping,
                         '-s', s_cluster_size])
    subprocess.check_call(map(str, seqnoise_cmd), cwd=cwd)
    return {"sequence_map": os.path.join(cwd, stub),
            "denoised_sequences": os.path.join(cwd, stub + '_cd.fa'),
            "mapping_file": os.path.join(cwd, stub + ".mapping")}

if __name__ == '__main__':
    main()
