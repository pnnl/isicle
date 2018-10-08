from os.path import *
from resources.nwchem.parseOutputs import XYZtoMFJ

# snakemake configuration
configfile: 'isicle/config.yaml'
include: 'molecular_dynamics.snakefile'


rule copyOver:
    input:
        join(config['path'], 'output', 'md', 'downselected', '{id}_{adduct}_{cycle}_{selected}.xyz')
    output:
        join(config['path'], 'output', 'dft', '{id}_{adduct}', 'cycle_{cycle}_{selected}', '{id}_{adduct}_{cycle}_{selected}.xyz')
    # group:
    #     'dft'
    shell:
        'cp {input} {output}'

# create .nw files based on template (resources/nwchem/template.nw)
rule createNW:
    input:
        rules.copyOver.output
    output:
        join(config['path'], 'output', 'dft', '{id}_{adduct}', 'cycle_{cycle}_{selected}', '{id}_{adduct}_{cycle}_{selected}.nw')
    # group:
    #     'dft'
    shell:
        'python isicle/resources/nwchem/generateNW.py {input} --template {config[nwchem][dft_template]}'

# run NWChem
rule NWChem:
    input:
        xyz = rules.createNW.input,
        nw = rules.createNW.output
    output:
        join(config['path'], 'output', 'dft', '{id}_{adduct}', 'cycle_{cycle}_{selected}', '{id}_{adduct}_{cycle}_{selected}.out')
    benchmark:
        join(config['path'], 'output', 'dft', 'benchmarks', '{id}_{adduct}_{cycle}_{selected}.dft.benchmark')
    # group:
    #     'dft'
    shell:
        '{config[nwchem][runscript]} {input.nw}'

# parse nwchem outputs (geometry files)
rule parseNWChem:
    input:
        rules.NWChem.output
    output:
        geom1 = join(config['path'], 'output', 'mobility', '{id}_{adduct}_{cycle}_{selected}_charge.mfj'),
        geom2 = join(config['path'], 'output', 'mobility', '{id}_{adduct}_{cycle}_{selected}_geom+charge.mfj'),
        charge1 = join(config['path'], 'output', 'mobility', '{id}_{adduct}_{cycle}_{selected}_charge.energy'),
        charge2 = join(config['path'], 'output', 'mobility', '{id}_{adduct}_{cycle}_{selected}_geom+charge.energy')
    benchmark:
        join(config['path'], 'output', 'nwchem', 'benchmarks', '{id}_{adduct}_{cycle}_{selected}.parse.benchmark')
    # group:
    #     'dft'
    run:
        XYZtoMFJ(input[0], join(config['path'], 'output', 'mobcal'))
