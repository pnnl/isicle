from os.path import *
from isicle.utils import cycles

# snakemake configuration
include: 'molecular_dynamics.snakefile'


rule copyOver:
    input:
        join(config['path'], 'output', 'md', 'downselected', '{id}_Ne_{cycle}_{selected}.xyz')
    output:
        join(config['path'], 'output', 'shielding', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.xyz')
    log:
        join(config['path'], 'output', 'shielding', 'logs', '{id}_{cycle}_{selected}.copy.log')
    benchmark:
        join(config['path'], 'output', 'shielding', 'benchmarks', '{id}_{cycle}_{selected}.copy.benchmark')
    # group:
    #     'shielding'
    shell:
        'cp {input} {output} &> {log}'


# create .nw files based on template (resources/nwchem/template.nw)
rule createShieldingConfig:
    input:
        rules.copyOver.output
    output:
        join(config['path'], 'output', 'shielding', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.nw')
    version:
        'python -m isicle.scripts.generateNW --version'
    log:
        join(config['path'], 'output', 'shielding', 'logs', '{id}_{cycle}_{selected}.create.log')
    benchmark:
        join(config['path'], 'output', 'shielding', 'benchmarks', '{id}_{cycle}_{selected}.create.benchmark')
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
        join(config['path'], 'output', 'shielding', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.out')
    version:
        "nwchem /dev/null | grep '(NWChem)' | awk '{print $6}'"
    log:
        join(config['path'], 'output', 'shielding', 'logs', '{id}_{cycle}_{selected}.nwchem.log')
    benchmark:
        join(config['path'], 'output', 'shielding', 'benchmarks', '{id}_{cycle}_{selected}.nwchem.benchmark')
    # group:
    #     'shielding'
    shell:
        'srun --mpi=pmi2 nwchem {input} > {output} 2> {log}'


# placeholder until shielding is done
rule parseShielding:
    input:
        rules.shielding.output
    output:
        shifts = join(config['path'], 'output', 'shifts', 'runs', '{id}_{cycle}_{selected}.shifts'),
        energy = join(config['path'], 'output', 'shifts', 'runs', '{id}_{cycle}_{selected}.energy')
    log:
        join(config['path'], 'output', 'shifts', 'runs', 'logs', '{id}_{cycle}_{selected}.log')
    benchmark:
        join(config['path'], 'output', 'shifts', 'runs', 'benchmarks', '{id}_{cycle}_{selected}.benchmark')
    # group:
    #     'shielding'
    run:
        outdir = dirname(output.shifts)
        shell('python -m isicle.scripts.parse_nwchem {input} %s --shielding &> {log}' % outdir)


rule combine:
    input:
        shifts = expand(join(config['path'], 'output', 'shifts', 'runs', '{{id}}_{cycle}_{selected}.shifts'),
                        cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2']),
        energy = expand(join(config['path'], 'output', 'shifts', 'runs', '{{id}}_{cycle}_{selected}.energy'),
                        cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2'])
    output:
        join(config['path'], 'output', 'shifts', 'conformer_shifts', '{id}.tsv')
    version:
        'python -m isicle.scripts.combine_shifts --version'
    log:
        join(config['path'], 'output', 'shifts', 'conformer_shifts', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'shifts', 'conformer_shifts', 'benchmarks', '{id}.benchmark')
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.combine_shifts {output} --infiles {input.shifts} --efiles {input.energy}  &> {log}'


# boltzmann averaging
rule boltzmannAverage:
    input:
        rules.combine.output
    output:
        join(config['path'], 'output', 'shifts', 'boltzmann_shifts', '{id}.tsv')
    version:
        'python -m isicle.scripts.boltzmann --version'
    log:
        join(config['path'], 'output', 'shifts', 'boltzmann_shifts', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'shifts', 'boltzmann_shifts', 'benchmarks', '{id}.benchmark')
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.boltzmann {input} {output} &> {log}'
