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
        expand(join(config['path'], 'output', 'pKa', '{id}.pka'), id=IDS)
        # expand(join(config['path'], 'output', 'pKa', '{id}_{adduct}.xyz'), id=IDS, adduct=ADDUCTS)

rule desalt:
    input:
        join(config['path'], 'input', '{id}.inchi')
    output:
        join(config['path'], 'output', 'desalted', '{id}.inchi')
    run:
        inchi = read_string(input[0])
        inchi = desalt(inchi)
        if inchi.startswith('InChI='):
            write_string(inchi, output[0])
        else:
            # raise IOError('OpenBabel', inchi)
            sys.exit(1)

rule neutralize:
    input:
        rules.desalt.output
    output:
        join(config['path'], 'output', 'neutralized', '{id}.inchi')
    run:
        inchi = read_string(input[0])
        inchi = neutralize(inchi)
        if inchi.startswith('InChI='):
            write_string(inchi, output[0])
        else:
            # raise IOError('neutralize', inchi)
            sys.exit(1)

rule majorTautomer:
    input:
        rules.neutralize.output
    output:
        join(config['path'], 'output', 'tautomer', '{id}.inchi')
    run:
        inchi = read_string(input[0])
        inchi = major_tautomer(inchi)
        if inchi.startswith('InChI='):
            write_string(inchi, output[0])
        else:
            # raise IOError('cxcalc', inchi)
            sys.exit(1)

rule calculateFormula:
    input:
        rules.majorTautomer.output
    output:
        join(config['path'], 'output', 'formula', '{id}.formula')
    run:
        inchi = read_string(input[0])
        formula = inchi2formula(inchi)
        write_string(formula, output[0])

rule calculateMass:
    input:
        rules.calculateFormula.output
    output:
        join(config['path'], 'output', 'mass', '{id}.mass')
    shell:
        'python resources/molmass.py `cat {input}` > {output}'

rule inchi2geom:
    input:
        inchi = rules.majorTautomer.output,
        formula = rules.calculateFormula.output,
        mass = rules.calculateMass.output
    output:
        mol = join(config['path'], 'output', 'parent_structures', 'mol', '{id}.mol'),
        png = join(config['path'], 'output', 'parent_structures', 'png', '{id}.png')
    run:
        inchi2geom(read_string(input[0]), output[0], output[1], ffield='gaff')

rule calculatepKa:
    input:
        rules.inchi2geom.output.mol
    output:
        join(config['path'], 'output', 'pKa', '{id}.pka')
    shell:
        'cxcalc pka -i -40 -x 40 -d large {input} > {output}'

# rule generateSodiated:
#     input:
#         geometry = rules.inchi2geom.output.mol,
#         pKa_values = rules.calculatepKa.output
#     output:
#         join(config['path'], 'output', 'adduct_structures', 'xyz', '{id}_+Na.xyz')
#         join(config['path'], 'output', 'adduct_structures', 'mol2', '{id}_+Na.mol2')
#     run:
#         # resources/adduct_creation/create_adducts.py

#         # Final output once we have CCS
#         # InChI, Processed InChI, Formula, Mass, CCS
