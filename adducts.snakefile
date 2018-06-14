from os.path import *
from resources.utils import *
import sys

# snakemake configuration
configfile: 'config.yaml'
localrules: all, desalt, neutralize, calculateMass

# infer wildcards from inputs, select adducts
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))
ADDUCTS = ['+Na', '+H', '-H']

# pseudo-rule that collects the target files
rule all:
    input:
        expand(join(config['path'], 'output', '3a_pKa', '{id}.pka'), id=IDS)
        # expand(join(config['path'], 'output', 'pKa', '{id}_{adduct}.xyz'), id=IDS, adduct=ADDUCTS)

# TODO: rule inchisToKeys:

rule desalt:
    input:
        join(config['path'], 'input', '{id}.inchi')
    output:
        join(config['path'], 'output', '0_desalted', '{id}.inchi')
    log:
        join(config['path'], 'output', '0_desalted', 'logs', '{id}.log')
    run:
        inchi = read_string(input[0])
        inchi = desalt(inchi, log[0])
        if inchi is not None:
            write_string(inchi, output[0])
        else:
            sys.exit(1)

rule neutralize:
    input:
        rules.desalt.output
    output:
        join(config['path'], 'output', '1_neutralized', '{id}.inchi')
    run:
        inchi = read_string(input[0])
        inchi = neutralize(inchi)
        if inchi is not None:
            write_string(inchi, output[0])
        else:
            sys.exit(1)

rule tautomerize:
    input:
        rules.neutralize.output
    output:
        join(config['path'], 'output', '2_tautomer', '{id}.inchi')
    log:
        join(config['path'], 'output', '2_tautomer', 'logs', '{id}.log')
    run:
        inchi = read_string(input[0])
        inchi = tautomerize(inchi, log=log[0])
        if inchi is not None:
            write_string(inchi, output[0])
        else:
            sys.exit(1)

rule calculateFormula:
    input:
        rules.tautomerize.output
    output:
        join(config['path'], 'output', '2a_formula', '{id}.formula')
    log:
        join(config['path'], 'output', '2a_formula', 'logs', '{id}.log')
    run:
        inchi = read_string(input[0])
        formula = inchi2formula(inchi, log=log[0])
        write_string(formula, output[0])

rule calculateMass:
    input:
        rules.calculateFormula.output
    output:
        join(config['path'], 'output', '2b_mass', '{id}.mass')
    shell:
        'python resources/molmass.py `cat {input}` > {output}'

rule inchi2geom:
    input:
        inchi = rules.tautomerize.output,
        formula = rules.calculateFormula.output,
        mass = rules.calculateMass.output
    output:
        mol = join(config['path'], 'output', '3_parent_structures', 'mol', '{id}.mol'),
        png = join(config['path'], 'output', '3_parent_structures', 'png', '{id}.png')
    run:
        inchi2geom(read_string(input[0]), output[0], output[1], ffield='gaff')

rule calculatepKa:
    input:
        rules.inchi2geom.output.mol
    output:
        join(config['path'], 'output', '3a_pKa', '{id}.pka')
    shell:
        'cxcalc pka -i -40 -x 40 -d large {input} > {output}'

# rule generateAdducts:
#     input:
#         geometry = rules.inchi2geom.output.mol,
#         pKa_values = rules.calculatepKa.output
#     output:
#         xyz = expand(join(config['path'], 'output', 'adduct_structures', 'xyz', {id}_{adduct}.xyz'), id=IDS, adduct=ADDUCTS),
#         mol2 = expand(join(config['path'], 'output', 'adduct_structures', 'mol2', {id}_{adduct}.mol2'), id=IDS, adduct=ADDUCTS)
#     run:
#         # resources/adduct_creation/create_adducts.py

#         # Final output once we have CCS
#         # InChI, Processed InChI, Formula, Mass, CCS
