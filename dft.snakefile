from os.path import *
from resources.nwchem.parseOutputs import XYZtoMFJ

# snakemake configuration
configfile: 'config.yaml'
include: 'molecular_dynamics.snakefile'
localrules: createNW, parseNWChem

# create .nw files based on template (resources/nwchem/template.nw)
rule createNW:
    input:
        rules.downselect.output.selected
    output:
        join(config['path'], 'output', 'nwchem', '{id}_{adduct}_{cycle}_{selected}.nw')
    group:
        'dft'
    shell:
        'python resources/nwchem/generateNW.py {input} --template {config[nwchem][template]}'

# run NWChem
rule NWChem:
    input:
        xyz = rules.createNW.input,
        nw = rules.createNW.output
    output:
        join(config['path'], 'output', 'nwchem', '{id}_{adduct}_{cycle}_{selected}.out')
    group:
        'dft'
    shell:
        '{config[nwchem][runscript]} {input.nw}'

# parse nwchem outputs (geometry files)
rule parseNWChem:
    input:
        rules.NWChem.output
    output:
        geom1 = join(config['path'], 'output', 'mobcal', '{id}_{adduct}_{cycle}_{selected}_charge.mfj'),
        geom2 = join(config['path'], 'output', 'mobcal', '{id}_{adduct}_{cycle}_{selected}_geom+charge.mfj'),
        charge1 = join(config['path'], 'output', 'mobcal', '{id}_{adduct}_{cycle}_{selected}_charge.energy'),
        charge2 = join(config['path'], 'output', 'mobcal', '{id}_{adduct}_{cycle}_{selected}_geom+charge.energy')
    group:
        'dft'
    run:
        XYZtoMFJ(input[0], join(config['path'], 'output', 'mobcal'))
