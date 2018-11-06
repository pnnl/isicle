from os.path import *
from pkg_resources import resource_filename


# snakemake configuration
include: 'dft.snakefile'
MOBCAL = resource_filename('isicle', 'resources/mobcal/%s' % config['mobcal']['exe'])

if config['mobcal']['params'] == 'default':
    PARAMS = resource_filename('isicle', 'resources/mobcal/mobcal.params')
else:
    PARAMS = config['mobcal']['params']

if config['mobcal']['atomtypes'] == 'default':
    ATOMS = resource_filename('isicle', 'resources/mobcal/atomtype_parameters.in')
else:
    ATOMS = config['mobcal']['atomtypes']


# run mobcal on geom+charge nwchem output
rule mobcal:
    input:
        rules.parseDFT.output.mfj
    output:
        join('output', 'mobility', 'mobcal', 'runs', '{id}_{adduct}_{cycle}_{selected}.out')
    version:
        '0.1.0'
    log:
        join('output', 'mobility', 'mobcal', 'runs', 'logs', '{id}_{adduct}_{cycle}_{selected}.log')
    benchmark:
        join('output', 'mobility', 'mobcal', 'runs', 'benchmarks', '{id}_{adduct}_{cycle}_{selected}.benchmark')
    # group:
    #     'mobility'
    shell:
        '{MOBCAL} {PARAMS} {ATOMS} {input} {output} &> {log}'


# parse mobcal output
rule parseMobcal:
    input:
        geom = expand(join('output', 'mobility', 'mobcal', 'runs', '{{id}}_{{adduct}}_{cycle}_{selected}.out'),
                      cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2']),
        energy = expand(join('output', 'mobility', 'mobcal', 'runs', '{{id}}_{{adduct}}_{cycle}_{selected}.energy'),
                        cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2'])
    output:
        join('output', 'mobility', 'mobcal', 'conformer_ccs', '{id}_{adduct}.tsv')
    version:
        'isicle --version'
    log:
        join('output', 'mobility', 'mobcal', 'conformer_ccs', 'logs', '{id}_{adduct}.log')
    benchmark:
        join('output', 'mobility', 'mobcal', 'conformer_ccs', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'mobility'
    shell:
        'python -m isicle.scripts.parse_mobcal {output} --infiles {input.geom} --efiles {input.energy} &> {log}'


# boltzmann averaging
rule boltzmannAverage:
    input:
        rules.parseMobcal.output
    output:
        join('output', 'mobility', 'mobcal', 'boltzmann_ccs', '{id}_{adduct}.tsv')
    version:
        'isicle --version'
    log:
        join('output', 'mobility', 'mobcal', 'boltzmann_ccs', 'logs', '{id}_{adduct}.log')
    benchmark:
        join('output', 'mobility', 'mobcal', 'boltzmann_ccs', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'mobility'
    shell:
        'python -m isicle.scripts.boltzmann {input} {output} --ccs &> {log}'


rule calibrate:
    input:
        rules.boltzmannAverage.output
    output:
        join('output', 'mobility', 'mobcal', 'calibrated_ccs', '{id}_{adduct}.tsv')
    version:
        'isicle --version'
    log:
        join('output', 'mobility', 'mobcal', 'calibrated_ccs', 'logs', '{id}_{adduct}.log')
    benchmark:
        join('output', 'mobility', 'mobcal', 'calibrated_ccs', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'mobility'
    run:
        m = config['ccs']['calibration'][wildcards.adduct]['m']
        b = config['ccs']['calibration'][wildcards.adduct]['b']
        shell('python -m isicle.scripts.calibrate {input} {output} {m} {b} &> {log}')
