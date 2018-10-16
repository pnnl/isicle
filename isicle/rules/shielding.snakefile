from os.path import *
from isicle.utils import cycles

# snakemake configuration
include: 'molecular_dynamics.snakefile'


rule copyOver:
    input:
        join(config['path'], 'output', 'md', 'downselected', '{id}_Ne_{cycle}_{selected}.xyz')
    output:
        join(config['path'], 'output', 'shielding', 'nwchem', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.xyz')
    log:
        join(config['path'], 'output', 'shielding', 'nwchem', 'logs', '{id}_{cycle}_{selected}.copy.log')
    benchmark:
        join(config['path'], 'output', 'shielding', 'nwchem', 'benchmarks', '{id}_{cycle}_{selected}.copy.benchmark')
    # group:
    #     'shielding'
    shell:
        'cp {input} {output} &> {log}'


# create .nw files based on template (resources/nwchem/template.nw)
rule createShieldingConfig:
    input:
        rules.copyOver.output
    output:
        join(config['path'], 'output', 'shielding', 'nwchem', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.nw')
    version:
        'python -m isicle.scripts.generateNW --version'
    log:
        join(config['path'], 'output', 'shielding', 'nwchem', 'logs', '{id}_{cycle}_{selected}.create.log')
    benchmark:
        join(config['path'], 'output', 'shielding', 'nwchem', 'benchmarks', '{id}_{cycle}_{selected}.create.benchmark')
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.generateNW {input} --shielding --shifts {config[nwchem][shifts]} \
         --template {config[nwchem][shielding_template]} &> {log}'


# run NWChem
rule shielding:
    input:
        rules.createShieldingConfig.output
    output:
        join(config['path'], 'output', 'shielding', 'nwchem', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.out')
    version:
        "nwchem /dev/null | grep '(NWChem)' | awk '{print $6}'"
    log:
        join(config['path'], 'output', 'shielding', 'nwchem', 'logs', '{id}_{cycle}_{selected}.nwchem.log')
    benchmark:
        join(config['path'], 'output', 'shielding', 'nwchem', 'benchmarks', '{id}_{cycle}_{selected}.nwchem.benchmark')
    # group:
    #     'shielding'
    shell:
        'srun --mpi=pmi2 nwchem {input} > {output} 2> {log}'


# placeholder until shielding is done
rule parseShielding:
    input:
        rules.shielding.output
    output:
        shielding = join(config['path'], 'output', 'shielding', 'parsed', '{id}_{cycle}_{selected}.shielding')
    log:
        join(config['path'], 'output', 'shielding', 'parsed', 'logs', '{id}_{cycle}_{selected}.log')
    benchmark:
        join(config['path'], 'output', 'shielding', 'parsed', 'benchmarks', '{id}_{cycle}_{selected}.benchmark')
    # group:
    #     'shielding'
    run:
        outdir = dirname(output.shielding)
        shell('python -m isicle.scripts.parse_nwchem {input} %s --shielding &> {log}' % outdir)


rule combine:
    input:
        expand(join(config['path'], 'output', 'shielding', 'parsed', '{{id}}_{cycle}_{selected}.shielding'),
               cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2'])
    output:
        join(config['path'], 'output', 'shielding', 'conformer_shielding', '{id}.tsv')
    version:
        'python -m isicle.scripts.combine_shifts --version'
    log:
        join(config['path'], 'output', 'shielding', 'conformer_shielding', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'shielding', 'conformer_shielding', 'benchmarks', '{id}.benchmark')
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.combine_shielding {input} {output} &> {log}'


# boltzmann averaging
rule boltzmannAverage:
    input:
        rules.combine.output
    output:
        join(config['path'], 'output', 'shielding', 'boltzmann_shielding', '{id}.tsv')
    version:
        'python -m isicle.scripts.boltzmann --version'
    log:
        join(config['path'], 'output', 'shielding', 'boltzmann_shielding', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'shielding', 'boltzmann_shielding', 'benchmarks', '{id}.benchmark')
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.boltzmann {input} {output} &> {log}'


rule shifts:
    input:
        shielding = rules.boltzmannAverage.output,
        ref = join(config['path'], 'output', 'shielding', 'boltzmann_shielding', '{config[nwchem][reference}.tsv')
    output:
        join(config['path'], 'output', 'shifts', '{id}.tsv')
    version:
        ''
    log:
        join(config['path'], 'output', 'shifts', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'shifts', 'benchmarks', '{id}.benchmark')
    shell:
        'echo'
