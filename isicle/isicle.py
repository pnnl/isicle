import argparse
from multiprocessing import cpu_count
from isicle import __version__
from isicle.utils import inchi2key, smi2key, write_string
from os.path import *
import pandas as pd
import os


def process(infile, outdir):
    df = pd.read_csv(infile, sep='\n', header=None)

    if not exists(outdir):
        os.makedirs(outdir)

    for row in df.values:
        mol = row[0]
        if 'InChI' in mol:
            key = inchi2key(mol)
            ext = '.inchi'
        else:
            key = smi2key(mol)
            ext = '.smi'

        if key is not None:
            write_string(mol, join(outdir, key + ext))
        else:
            print('Failed to hash %s.' % mol)


def cli():
    p = {}
    # overall
    p['global'] = argparse.ArgumentParser(description='ISiCLE: in silico chemical library engine')
    p['global'].add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    # parent with shared args
    p['parent'] = argparse.ArgumentParser(add_help=False)
    p['parent'].add_argument('--config', metavar='PATH', default='config.yaml', help='path to isicle yaml configuration file')

    # parent with snakemake args
    p['sm'] = argparse.ArgumentParser(add_help=False)
    p['sm'].add_argument('--cluster', metavar='PATH', help='path to cluster execution yaml configuration file')
    p['sm'].add_argument('--dryrun', action='store_true', help='perform a dry run')
    p['sm'].add_argument('--unlock', action='store_true', help='unlock directory')
    p['sm'].add_argument('--cores', metavar='N', type=int, default=cpu_count(),
                         help='number of cores used for execution (local execution only)')
    p['sm'].add_argument('--jobs', metavar='N', type=int, default=1000,
                         help='number of simultaneous jobs to submit to a slurm queue (cluster execution only)')

    # subparsers
    p['subparsers'] = p['global'].add_subparsers(title='commands', dest='which')

    # prep mode
    p['prep'] = p['subparsers'].add_parser('prep', parents=[p['parent']],
                                           description='ISiCLE input preparation module',
                                           help='input preparation module')
    p['prep'].add_argument('infile', help='path to input file containing inchi and/or smiles strings')

    # ccs mode
    p['ccs'] = p['subparsers'].add_parser('ccs', parents=[p['parent'], p['sm']],
                                          description='ISiCLE CCS calculation module',
                                          help='ccs calculation module')
    p['ccs'].add_argument('mode', choices=['standard', 'lite'], help='ccs calculation mode')

    # shifts mode
    p['shifts'] = p['subparsers'].add_parser('shifts', parents=[p['parent'], p['sm']],
                                             description='ISiCLE NMR chemical shifts calculation module',
                                             help='chemical shift calculation module')

    args = p['global'].parse_args()

    import subprocess
    from pkg_resources import resource_filename

    # check for config
    if not isfile(args.config):
        p['global'].error('Snakemake YAML configuration file not found.')

    # input processing
    if args.which == 'prep':
        import yaml

        with open(args.config, 'r') as f:
            outdir = join(yaml.load(f)['path'], 'input')

        process(args.infile, outdir)

    else:
        # ccs module
        if args.which == 'ccs':
            # standard mode
            if args.mode == 'standard':
                cmd = 'snakemake --snakefile %s --configfile %s -k --rerun-incomplete' % \
                      (resource_filename('isicle', 'rules/ccs_standard.snakefile'), args.config)
            # lite mode
            elif args.mode == 'lite':
                cmd = 'snakemake --snakefile %s --configfile %s -k --rerun-incomplete' % \
                      (resource_filename('isicle', 'rules/ccs_lite.snakefile'), args.config)

        # chemical shifts module
        elif args.which == 'shifts':
            cmd = 'snakemake --snakefile %s --configfile %s -k --rerun-incomplete' % \
                  (resource_filename('isicle', 'rules/shifts.snakefile'), args.config)

        # cluster configuration
        if args.cluster is not None:
            cmd += ' --cluster-config %s' % args.cluster
            cmd += ' --cluster "sbatch -A {cluster.account} -N {cluster.nodes} -t {cluster.time} -J {cluster.name} --ntasks-per-node {cluster.ntasks}"'
            cmd += ' -j %s' % args.jobs
        # local execution
        else:
            cmd += ' --cores %s' % args.cores

        # additional options
        if args.dryrun:
            cmd += ' --dryrun'
        if args.unlock:
            cmd += ' --unlock'

        # execute
        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    cli()
