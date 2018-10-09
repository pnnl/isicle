from os.path import *

# snakemake configuration
include: 'dft.snakefile'


# run mobcal on geom+charge nwchem output
rule mobcal:
    input:
        rules.parseNWChem.output.geom2
    output:
        join(config['path'], 'output', 'mobility', 'mobcal', 'runs', '{id}_{adduct}_{cycle}_{selected}_geom+charge.out')
    benchmark:
        join(config['path'], 'output', 'mobility', 'mobcal', 'runs', 'benchmarks', '{id}_{adduct}_{cycle}_{selected}.benchmark')
    # group:
    #     'mobility'
    shell:
        '{config[mobcal][exe]} {config[mobcal][params]} {config[mobcal][atomtypes]} {input} {output}'


# parse mobcal output
rule parseMobcal:
    input:
        geom = expand(join(config['path'], 'output', 'mobility', 'mobcal', 'runs', '{{id}}_{{adduct}}_{cycle}_{selected}_geom+charge.out'),
                      cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2']),
        energy = expand(join(config['path'], 'output', 'mobility', 'mobcal', 'runs', '{{id}}_{{adduct}}_{cycle}_{selected}_geom+charge.energy'),
                        cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2'])
    output:
        join(config['path'], 'output', 'mobility', 'mobcal', 'conformer_ccs', '{id}_{adduct}.tsv')
    log:
        join(config['path'], 'output', 'mobility', 'mobcal', 'conformer_ccs', 'logs', '{id}_{adduct}.log')
    benchmark:
        join(config['path'], 'output', 'mobility', 'mobcal', 'conformer_ccs', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'mobility'
    shell:
        'python isicle/parse_mobcal.py {input.geom} {input.energy} {output} &> {log}'


# boltzmann averaging
rule boltzmannAverage:
    input:
        rules.parseMobcal.output
    output:
        join(config['path'], 'output', 'mobility', 'mobcal', 'boltzmann_ccs', '{id}_{adduct}.tsv')
    log:
        join(config['path'], 'output', 'mobility', 'mobcal', 'boltzmann_ccs', 'logs', '{id}_{adduct}.log')
    benchmark:
        join(config['path'], 'output', 'mobility', 'mobcal', 'boltzmann_ccs', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'mobility'
    shell:
        'python isicle/boltzmann.py {input} {output} &> {log}'
