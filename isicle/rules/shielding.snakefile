from os.path import *
import json
from isicle.utils import cycles

# snakemake configuration
include: 'molecular_dynamics.snakefile'


rule copyOver:
    input:
        abspath(join('output', 'md', 'downselected', '{id}_neutral_{cycle}_{selected}.xyz'))
    output:
        abspath(join('output', 'shielding', 'nwchem', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.xyz'))
    log:
        abspath(join('output', 'shielding', 'nwchem', 'logs', '{id}_{cycle}_{selected}.copy.log'))
    benchmark:
        abspath(join('output', 'shielding', 'nwchem', 'benchmarks', '{id}_{cycle}_{selected}.copy.benchmark'))
    # group:
    #     'shielding'
    shell:
        'cp {input} {output} &> {log}'


# create .nw files based on template (resources/nwchem/template.nw)
rule createShieldingConfig:
    input:
        rules.copyOver.output
    output:
        abspath(join('output', 'shielding', 'nwchem', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.nw'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'shielding', 'nwchem', 'logs', '{id}_{cycle}_{selected}.create.log'))
    benchmark:
        abspath(join('output', 'shielding', 'nwchem', 'benchmarks', '{id}_{cycle}_{selected}.create.benchmark'))
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.generateNW {input} --shielding --shifts {config[nwchem][shifts]} \
         --template {config[nwchem][shielding_template]} --solvent {config[nwchem][solvent]} &> {log}'


# run NWChem
rule shielding:
    input:
        rules.createShieldingConfig.output
    output:
        abspath(join('output', 'shielding', 'nwchem', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.out'))
    version:
        "nwchem /dev/null | grep '(NWChem)' | awk '{print $6}'"
    log:
        abspath(join('output', 'shielding', 'nwchem', 'logs', '{id}_{cycle}_{selected}.nwchem.log'))
    benchmark:
        abspath(join('output', 'shielding', 'nwchem', 'benchmarks', '{id}_{cycle}_{selected}.nwchem.benchmark'))
    # group:
    #     'shielding'
    shell:
        'nwchem {input} > {output} 2> {log}'


# placeholder until shielding is done
rule parseShielding:
    input:
        rules.shielding.output
    output:
        shielding = abspath(join('output', 'shielding', 'parsed', '{id}_{cycle}_{selected}.shielding'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'shielding', 'parsed', 'logs', '{id}_{cycle}_{selected}.log'))
    benchmark:
        abspath(join('output', 'shielding', 'parsed', 'benchmarks', '{id}_{cycle}_{selected}.benchmark'))
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.parse_nwchem {input} {output.shielding} --shielding &> {log}'


rule combine:
    input:
        expand(abspath(join('output', 'shielding', 'parsed', '{{id}}_{cycle}_{selected}.shielding')),
               cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2'])
    output:
        abspath(join('output', 'shielding', 'conformer_shielding', '{id}.tsv'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'shielding', 'conformer_shielding', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'shielding', 'conformer_shielding', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.combine_shielding {input} {output} &> {log}'


# boltzmann averaging
rule boltzmannAverage:
    input:
        rules.combine.output
    output:
        abspath(join('output', 'shielding', 'boltzmann_shielding', '{id}.tsv'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'shielding', 'boltzmann_shielding', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'shielding', 'boltzmann_shielding', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'shielding'
    shell:
        'python -m isicle.scripts.boltzmann {input} {output} --shielding &> {log}'


rule shifts:
    input:
        shielding = rules.boltzmannAverage.output
    output:
        abspath(join('output', 'shifts', '{id}.tsv'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'shifts', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'shifts', 'benchmarks', '{id}.benchmark'))
    run:
        ref = json.dumps(config['nwchem']['reference'])
        shell("python -m isicle.scripts.calculate_shifts {input.shielding} '{ref}' {output} &> {log}")
