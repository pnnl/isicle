from os.path import *
from isicle.utils import cycles

# snakemake configuration
include: 'molecular_dynamics.snakefile'


rule copyOver:
    input:
        join('output', 'md', 'downselected', '{id}_neutral_{cycle}_{selected}.xyz')
    output:
        join('output', 'shielding', 'nwchem', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.xyz')
    log:
        join('output', 'shielding', 'nwchem', 'logs', '{id}_{cycle}_{selected}.copy.log')
    benchmark:
        join('output', 'shielding', 'nwchem', 'benchmarks', '{id}_{cycle}_{selected}.copy.benchmark')
    # group:
    #     'shielding'
    shell:
        'cp {input} {output} &> {log}'


# create .nw files based on template (resources/nwchem/template.nw)
rule createShieldingConfig:
    input:
        rules.copyOver.output
    output:
        join('output', 'shielding', 'nwchem', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.nw')
    version:
        'python -m isicle.scripts.generateNW --version'
    log:
        join('output', 'shielding', 'nwchem', 'logs', '{id}_{cycle}_{selected}.create.log')
    benchmark:
        join('output', 'shielding', 'nwchem', 'benchmarks', '{id}_{cycle}_{selected}.create.benchmark')
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
        join('output', 'shielding', 'nwchem', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.out')
    version:
        "nwchem /dev/null | grep '(NWChem)' | awk '{print $6}'"
    log:
        join('output', 'shielding', 'nwchem', 'logs', '{id}_{cycle}_{selected}.nwchem.log')
    benchmark:
        join('output', 'shielding', 'nwchem', 'benchmarks', '{id}_{cycle}_{selected}.nwchem.benchmark')
    # group:
    #     'shielding'
    run:
        shell('cd /scratch')
        shell('srun --mpi=pmi2 nwchem {input} > {output} 2> {log}')


# placeholder until shielding is done
rule parseShielding:
    input:
        rules.shielding.output
    output:
        shielding = join('output', 'shielding', 'parsed', '{id}_{cycle}_{selected}.shielding')
    log:
        join('output', 'shielding', 'parsed', 'logs', '{id}_{cycle}_{selected}.log')
    benchmark:
        join('output', 'shielding', 'parsed', 'benchmarks', '{id}_{cycle}_{selected}.benchmark')
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.parse_nwchem {input} {output.shielding} --shielding &> {log}'


rule combine:
    input:
        expand(join('output', 'shielding', 'parsed', '{{id}}_{cycle}_{selected}.shielding'),
               cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2'])
    output:
        join('output', 'shielding', 'conformer_shielding', '{id}.tsv')
    version:
        'python -m isicle.scripts.combine_shifts --version'
    log:
        join('output', 'shielding', 'conformer_shielding', 'logs', '{id}.log')
    benchmark:
        join('output', 'shielding', 'conformer_shielding', 'benchmarks', '{id}.benchmark')
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.combine_shielding {input} {output} &> {log}'


# boltzmann averaging
rule boltzmannAverage:
    input:
        rules.combine.output
    output:
        join('output', 'shielding', 'boltzmann_shielding', '{id}.tsv')
    version:
        'python -m isicle.scripts.boltzmann --version'
    log:
        join('output', 'shielding', 'boltzmann_shielding', 'logs', '{id}.log')
    benchmark:
        join('output', 'shielding', 'boltzmann_shielding', 'benchmarks', '{id}.benchmark')
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.boltzmann {input} {output} --shielding &> {log}'


rule shifts:
    input:
        shielding = rules.boltzmannAverage.output,
        ref = join('output', 'shielding', 'boltzmann_shielding', '%s.tsv' % config['nwchem']['reference'])
    output:
        join('output', 'shifts', '{id}.tsv')
    version:
        ''
    log:
        join('output', 'shifts', 'logs', '{id}.log')
    benchmark:
        join('output', 'shifts', 'benchmarks', '{id}.benchmark')
    shell:
        'python -m isicle.scripts.calculate_shifts {input.shielding} {input.ref} {output} &> {log}'
