import argparse
from multiprocessing import cpu_count
from os.path import *
import pandas as pd
import os
from pkg_resources import resource_filename
from snakemake import snakemake
from isicle.utils import inchi2key, smi2key, write_string
from isicle import __version__
from isicle import export
from isicle import preprocess


def process(infile):
    df = pd.read_csv(infile, sep='\n', header=None)
    df.dropna(inplace=True)

    if not exists('input'):
        os.mkdir('input')

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

    # parent with snakemake args
    p['sm'] = argparse.ArgumentParser(add_help=False)
    p['sm'].add_argument('--config', metavar='PATH', default='config.yaml', help='path to isicle yaml configuration file')
    p['sm'].add_argument('--dryrun', action='store_true', help='perform a dry run')
    p['sm'].add_argument('--unlock', action='store_true', help='unlock directory')
    p['sm'].add_argument('--touch', action='store_true', help='touch output files only')
    p['sm'].add_argument('--latency', metavar='N', type=int, default=3, help='specify filesystem latency (seconds)')
    p['sm'].add_argument('--cores', metavar='N', type=int, default=cpu_count(),
                         help='number of cores used for execution (local execution only)')
    p['sm'].add_argument('--count', metavar='N', type=int,
                         help='number of molecules to run (limits DAG size)')
    p['sm'].add_argument('--start', metavar='IDX', type=int, default=0,
                         help='starting molecule index (for use with --count)')

    # cluster-specific options
    p['clust'] = p['sm'].add_argument_group('cluster arguments')
    p['clust'].add_argument('--cluster', metavar='PATH', help='path to cluster execution yaml configuration file')
    p['clust'].add_argument('--jobs', metavar='N', type=int, default=1000,
                            help='number of simultaneous jobs to submit to a slurm queue')

    # subparsers
    p['subparsers'] = p['global'].add_subparsers(title='commands', dest='which')

    # prep mode
    p['prep'] = p['subparsers'].add_parser('prep',
                                           description='ISiCLE input preparation module',
                                           help='input preparation module')
    p['prep'].add_argument('infile', help='path to input file containing inchi and/or smiles strings')

    # ccs mode
    p['ccs'] = p['subparsers'].add_parser('ccs', parents=[p['sm']],
                                          description='ISiCLE CCS calculation module',
                                          help='ccs calculation module')
    p['ccs'].add_argument('mode', choices=['standard', 'lite'], help='ccs calculation mode')

    # shifts mode
    p['shifts'] = p['subparsers'].add_parser('shifts', parents=[p['sm']],
                                             description='ISiCLE NMR chemical shifts calculation module',
                                             help='chemical shift calculation module')

    # free energy mode
    p['energy'] = p['subparsers'].add_parser('energy', parents=[p['sm']],
                                             description="ISiCLE Gibb's Free Energy calculation module",
                                             help="gibb's free energy alculation module")

    # export
    p['export'] = p['subparsers'].add_parser('export',
                                             description="ISiCLE export module",
                                             help="result export module")

    p['process'] = p['subparsers'].add_parser('process',
                                              description="ISiCLE input processing module",
                                              help="input processing module")
    p['process'].add_argument('infile', help='path to input .tsv file containing smiles strings')
    p['process'].add_argument('outfile', help='path to output .tsv file')

    # export subparser
    p['export_subp'] = p['export'].add_subparsers(title='commands', dest='exportmode')

    # ccs
    p['export_ccs'] = p['export_subp'].add_parser('ccs',
                                                  description='Export CCS',
                                                  help='export ccs')
    p['export_ccs'].add_argument('mode', choices=['standard', 'lite'], help='ccs calculation mode')
    p['export_ccs'].add_argument('path', help='path to .tsv file for exported results')

    # shifts
    p['export_shifts'] = p['export_subp'].add_parser('shifts',
                                                     description='Export NMR chemical shifts',
                                                     help='export nmr chemical shifts')
    p['export_shifts'].add_argument('path', help='path to folder for exported results')

    # energy
    p['export_energy'] = p['export_subp'].add_parser('energy',
                                                     description="Export Gibb's free energy values",
                                                     help='export free energy values')
    p['export_energy'].add_argument('path', help='path to .tsv for exported results')

    args = p['global'].parse_args()

    # input processing
    if args.which == 'prep':
        process(args.infile)

    elif args.which == 'export':
        # ccs
        if args.exportmode == 'ccs':
            export.ccs(args.path, mode=args.mode)

        # shifts
        elif args.exportmode == 'shifts':
            export.shifts(args.path)

        # free energy
        elif args.exportmode == 'energy':
            export.energy(args.path)

        # none selected
        else:
            p['export'].print_help()

    elif args.which == 'process':
        df = pd.read_csv(args.infile, sep='\t')
        df['Processed'] = preprocess.process(df['SMILES'].values)
        df.to_csv(args.outfile, sep='\t', index=False)

    # simulation modules
    elif args.which in ['ccs', 'shifts', 'energy']:
        # check for config
        if not isfile(args.config):
            p['global'].error('Snakemake YAML configuration file not found.')

        # cluster config
        if args.cluster is not None:
            cluster = "sbatch -A {cluster.account} -N {cluster.nodes} -t {cluster.time} -J {cluster.name} --ntasks-per-node {cluster.ntasks}"
        else:
            cluster = None

        # start/stop config
        if args.count is not None:
            config = {'start': args.start, 'stop': args.start + args.count}
        else:
            config = {}

        # ccs module
        if args.which == 'ccs':
            # standard mode
            if args.mode == 'standard':
                snakemake(resource_filename('isicle', 'rules/ccs_standard.snakefile'),
                          configfiles=[args.config],
                          config=config,
                          cluster_config=args.cluster,
                          cluster=cluster,
                          keepgoing=True,
                          force_incomplete=True,
                          cores=args.cores,
                          nodes=args.jobs,
                          dryrun=args.dryrun,
                          unlock=args.unlock,
                          touch=args.touch,
                          latency_wait=args.latency)
            # lite mode
            elif args.mode == 'lite':
                snakemake(resource_filename('isicle', 'rules/ccs_lite.snakefile'),
                          configfiles=[args.config],
                          config=config,
                          cluster_config=args.cluster,
                          cluster=cluster,
                          keepgoing=True,
                          force_incomplete=True,
                          cores=args.cores,
                          nodes=args.jobs,
                          dryrun=args.dryrun,
                          unlock=args.unlock,
                          touch=args.touch,
                          latency_wait=args.latency)

        # chemical shifts module
        elif args.which == 'shifts':
            snakemake(resource_filename('isicle', 'rules/shifts.snakefile'),
                      configfiles=[args.config],
                      config=config,
                      cluster_config=args.cluster,
                      cluster=cluster,
                      keepgoing=True,
                      force_incomplete=True,
                      cores=args.cores,
                      nodes=args.jobs,
                      dryrun=args.dryrun,
                      unlock=args.unlock,
                      touch=args.touch,
                      latency_wait=args.latency)

        # gfe module
        elif args.which == 'energy':
            snakemake(resource_filename('isicle', 'rules/free_energy.snakefile'),
                      configfiles=[args.config],
                      config=config,
                      cluster_config=args.cluster,
                      cluster=cluster,
                      keepgoing=True,
                      force_incomplete=True,
                      cores=args.cores,
                      nodes=args.jobs,
                      dryrun=args.dryrun,
                      unlock=args.unlock,
                      touch=args.touch,
                      latency_wait=args.latency)

    # no module selected
    else:
        p['global'].print_help()


if __name__ == '__main__':
    cli()
