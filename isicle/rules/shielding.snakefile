from os.path import *
from helpers import cycles

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
    log:
        join(config['path'], 'output', 'shielding', 'logs', '{id}_{cycle}_{selected}.create.log')
    benchmark:
        join(config['path'], 'output', 'shielding', 'benchmarks', '{id}_{cycle}_{selected}.create.benchmark')
    # group:
    #     'shielding'
    shell:
        'python isicle/generateNW.py {input} --template {config[nwchem][shielding_template]} &> {log}'


# run NWChem
rule shielding:
    input:
        rules.createShieldingConfig.output
    output:
        join(config['path'], 'output', 'shielding', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.out')
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
        expand(join(config['path'], 'output', 'shielding', '{{id}}', 'cycle_{cycle}_{selected}', '{{id}}_{cycle}_{selected}.out'),
               selected=['s', 'd1', 'd2'], cycle=cycles(config['amber']['cycles']))
    output:
        join(config['path'], 'output', 'chemical_shifts', '{id}.tsv')
    log:
        join(config['path'], 'output', 'chemical_shifts', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'chemical_shifts', 'benchmarks', '{id}.benchmark')
    # group:
    #     'dft'
    # run:
    #     outdir = dirname(output.geom2)
    #     shell('python isicle/parse_nwchem.py {input} %s --mode shielding &> {log}' % outdir)
    shell:
        'touch {output}'
