import argparse
from multiprocessing import cpu_count
from isicle import __version__
from isicle.utils import inchi2key, smi2key, write_string
from os.path import *
import pandas as pd


def process(infile):
    df = pd.read_csv(infile, sep='\n', header=None)

    if not exists('input'):
        mkdir('input')

    for row in df.values:
        mol = row[0]
        if 'InChI' in mol:
            key = inchi2key(mol)
            ext = '.inchi'
        else:
            key = smi2key(mol)
            ext = '.smi'

        if key is not None:
            write_string(mol, join('input', key + ext))
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
    p['sm'].add_argument('--start', metavar='N', type=int, default=0,
                         help='starting molecule index')
    p['sm'].add_argument('--count', metavar='N', type=int,
                         help='number of molecules to run (limits DAG size)')

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

    from pkg_resources import resource_filename
    from snakemake import snakemake

    # check for config
    if not isfile(args.config):
        p['global'].error('Snakemake YAML configuration file not found.')

    # input processing
    if args.which == 'prep':
        process(args.infile)

    # simulation modules
    else:
        # cluster config
        if args.cluster is not None:
            cluster = "sbatch -A {cluster.account} -N {cluster.nodes} -t {cluster.time} -J {cluster.name} --ntasks-per-node {cluster.ntasks}"
        else:
            cluster = None

        # start/stop config
        if args.count is not None:
            config = {'start': args.start, 'stop': args.start + args.count}
        else:
            config = {'start': args.start, 'stop': None}

        # ccs module
        if args.which == 'ccs':
            # standard mode
            if args.mode == 'standard':
                snakemake(resource_filename('isicle', 'rules/ccs_standard.snakefile'),
                          configfile=args.config,
                          config=config,
                          cluster_config=args.cluster,
                          cluster=cluster,
                          keepgoing=True,
                          force_incomplete=True,
                          cores=args.cores,
                          nodes=args.jobs,
                          dryrun=args.dryrun,
                          unlock=args.unlock)
            # lite mode
            elif args.mode == 'lite':
                snakemake(resource_filename('isicle', 'rules/ccs_lite.snakefile'),
                          configfile=args.config,
                          config=config,
                          cluster_config=args.cluster,
                          cluster=cluster,
                          keepgoing=True,
                          force_incomplete=True,
                          cores=args.cores,
                          nodes=args.jobs,
                          dryrun=args.dryrun,
                          unlock=args.unlock)

        # chemical shifts module
        elif args.which == 'shifts':
            snakemake(resource_filename('isicle', 'rules/shifts.snakefile'),
                      configfile=args.config,
                      config=config,
                      cluster_config=args.cluster,
                      cluster=cluster,
                      keepgoing=True,
                      force_incomplete=True,
                      cores=args.cores,
                      nodes=args.jobs,
                      dryrun=args.dryrun,
                      unlock=args.unlock)


if __name__ == '__main__':
    cli()
