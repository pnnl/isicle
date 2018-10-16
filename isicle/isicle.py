import argparse
from multiprocessing import cpu_count
from isicle import __version__
from os.path import *


def cli():
    parser = argparse.ArgumentParser(description='Execute ISiCLE simulations.')
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    config = parser.add_argument_group('Snakemake configuration')
    config.add_argument('--config', required=True, help='Path to ISiCLE configuration file.')
    config.add_argument('--cluster-config', help='Path to cluster execution configuration file.')
    config.add_argument('--dryrun', action='store_true', help='Perform a dry run.')
    config.add_argument('--unlock', action='store_true', help='Unlock directory.')

    parallel = config.add_mutually_exclusive_group()
    parallel.add_argument('--cores', type=int, default=cpu_count(), help='Number of cores used for execution (local execution only).')
    parallel.add_argument('--jobs', type=int, default=1000, help='Number of simultaneous jobs to submit to a slurm queue (cluster execution only).')

    prop = parser.add_argument_group('Property calculation')
    mprop = prop.add_mutually_exclusive_group(required=True)
    mprop.add_argument('--ccs', action='store_true', help='Calculate collision cross section.')
    mprop.add_argument('--shifts', action='store_true', help='Calculate NMR chemical shifts.')

    mode = parser.add_argument_group('Calculation mode (CCS only)')
    mmode = mode.add_mutually_exclusive_group()
    mmode.add_argument('--standard', action='store_true', help='Standard calculation mode.')
    mmode.add_argument('--lite', action='store_true', help='Lite calculation mode.')

    args = parser.parse_args()

    import subprocess
    from pkg_resources import resource_filename

    if args.ccs is True:
        if args.standard is True:
            cmd = 'snakemake --snakefile %s --configfile %s -k --rerun-incomplete' % \
                  (resource_filename('isicle', 'rules/ccs_standard.snakefile'), args.config)
        elif args.lite is True:
            cmd = 'snakemake --snakefile %s --configfile %s -k --rerun-incomplete' % \
                  (resource_filename('isicle', 'rules/ccs_lite.snakefile'), args.config)
        else:
            parser.error('Please select a CCS calculation mode.')
    elif args.shifts is True:
        cmd = 'snakemake --snakefile %s --configfile %s -k --rerun-incomplete' % \
              (resource_filename('isicle', 'rules/shifts.snakefile'), args.config)

    if args.cluster_config is not None:
        cmd += ' --cluster-config %s' % args.cluster_config
        cmd += ' --cluster "sbatch -A {cluster.account} -N {cluster.nodes} -t {cluster.time} -J {cluster.name} --ntasks-per-node {cluster.ntasks}"'
        cmd += ' -j %s' % args.jobs
    else:
        cmd += ' --cores %s' % args.cores

    if args.dryrun:
        cmd += ' --dryrun'

    if args.unlock:
        cmd += ' --unlock'

    subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    cli()
