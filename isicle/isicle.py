import argparse
from multiprocessing import cpu_count
from isicle import __version__
from isicle.utils import inchi2key, smi2key
from os.path import *
import pandas as pd


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
    p['global'] = argparse.ArgumentParser(description='execute isicle simulations')
    p['subparsers'] = p['global'].add_subparsers(title='commands')
    p['global'].add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    p['parent'] = argparse.ArgumentParser(add_help=False)
    p['parent'].add_argument('config', help='path to isicle configuration file')

    p['sm'] = argparse.ArgumentParser(add_help=False)
    p['smc'] = p['sm'].add_argument_group('snakemake configuration')
    p['smc'].add_argument('--cluster', metavar='PATH', help='path to cluster execution configuration file')
    p['smc'].add_argument('--dryrun', action='store_true', help='perform a dry run')
    p['smc'].add_argument('--unlock', action='store_true', help='unlock directory')

    p['smp'] = p['smc'].add_mutually_exclusive_group()
    p['smp'].add_argument('--cores', metavar='N', type=int, default=cpu_count(), help='number of cores used for execution (local execution only)')
    p['smp'].add_argument('--jobs', metavar='N', type=int, default=1000, help='number of simultaneous jobs to submit to a slurm queue (cluster execution only)')

    p['prep'] = p['subparsers'].add_parser('prep', parents=[p['parent']], help='input preparation module')
    p['ccs'] = p['subparsers'].add_parser('ccs', parents=[p['parent'], p['sm']], help='ccs calculation module')
    p['shifts'] = p['subparsers'].add_parser('shifts', parents=[p['parent'], p['sm']], help='nmr chemical shifts calculation module')

    p['global'].parse_args()

    # # reorder default groups
    # optional = parser._action_groups.pop()
    # required = parser.add_argument_group('required arguments')
    # parser._action_groups.append(optional)

    # required.add_argument('mode', nargs='+', help='isicle operation mode (prep, ccs-standard, ccs-lite, shifts)')
    # required.add_argument('--config', metavar='PATH', required=True, help='path to isicle configuration file')
    # optional.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    # config = parser.add_argument_group('snakemake configuration')
    # config.add_argument('--cluster', metavar='PATH', help='path to cluster execution configuration file')
    # config.add_argument('--dryrun', action='store_true', help='perform a dry run')
    # config.add_argument('--unlock', action='store_true', help='unlock directory')

    # parallel = config.add_mutually_exclusive_group()
    # parallel.add_argument('--cores', metavar='N', type=int, default=cpu_count(), help='number of cores used for execution (local execution only)')
    # parallel.add_argument('--jobs', metavar='N', type=int, default=1000, help='number of simultaneous jobs to submit to a slurm queue (cluster execution only)')

    # args = parser.parse_args()

    # import subprocess
    # from pkg_resources import resource_filename

    # # input processing
    # if args.mode[0] == 'prep':
    #     import yaml

    #     with open(args.config, 'r') as f:
    #         outdir = join(yaml.load(f)['path'], 'input')

    #     process(args.mode[1], outdir)

    # else:
    #     # standard mode
    #     if args.mode[0] == 'ccs-standard':
    #         cmd = 'snakemake --snakefile %s --configfile %s -k --rerun-incomplete' % \
    #               (resource_filename('isicle', 'rules/ccs_standard.snakefile'), args.config)
    #     # lite mode
    #     elif args.mode[0] == 'ccs-lite':
    #         cmd = 'snakemake --snakefile %s --configfile %s -k --rerun-incomplete' % \
    #               (resource_filename('isicle', 'rules/ccs_lite.snakefile'), args.config)

    #     # chemical shifts module
    #     elif args.mode[0] == 'shifts':
    #         cmd = 'snakemake --snakefile %s --configfile %s -k --rerun-incomplete' % \
    #               (resource_filename('isicle', 'rules/shifts.snakefile'), args.config)

    #     # cluster configuration
    #     if args.cluster is not None:
    #         cmd += ' --cluster-config %s' % args.cluster
    #         cmd += ' --cluster "sbatch -A {cluster.account} -N {cluster.nodes} -t {cluster.time} -J {cluster.name} --ntasks-per-node {cluster.ntasks}"'
    #         cmd += ' -j %s' % args.jobs
    #     # local execution
    #     else:
    #         cmd += ' --cores %s' % args.cores

    #     # additional options
    #     if args.dryrun:
    #         cmd += ' --dryrun'
    #     if args.unlock:
    #         cmd += ' --unlock'

    #     # execute
    #     subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    cli()
