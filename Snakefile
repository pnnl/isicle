from os.path import *
from resources.utils import *
import sys

# snakemake configuration
configfile: 'config.yaml'
localrules: all, desalt, neutralize, calculateMass


# infer wildcards from inputs
IDS, = glob_wildcards(join(config['path'], 'input', 'molid_{id}.inchi'))

# pseudo-rule that collects the target files
rule all:
    input:
        expand(join(config['path'], 'output', 'pKa', 'molid_{id}_3D.pka'), id=IDS)

rule desalt:
    input:
        join(config['path'], 'input', 'molid_{id}.inchi')
    output:
        join(config['path'], 'output', 'desalted', 'molid_{id}.inchi')
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
        join(config['path'], 'output', 'neutralized', 'molid_{id}.inchi')
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
        join(config['path'], 'output', 'tautomer', 'molid_{id}.inchi')
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
        join(config['path'], 'output', 'formula', 'molid_{id}.formula')
    run:
        inchi = read_string(input[0])
        formula = inchi2formula(inchi)
        write_string(formula, output[0])

rule calculateMass:
    input:
        rules.calculateFormula.output
    output:
        join(config['path'], 'output', 'mass', 'molid_{id}.mass')
    shell:
        'python resources/molmass.py `cat {input}` > {output}'

rule inchi2geom:
    input:
        inchi = rules.majorTautomer.output,
        formula = rules.calculateFormula.output,
        mass = rules.calculateMass.output
    output:
        mol2D = join(config['path'], 'output', 'mol', 'molid_{id}_2D.mol'),
        mol3D = join(config['path'], 'output', 'mol', 'molid_{id}_3D.mol'),
        xyz = join(config['path'], 'output', 'xyz', 'molid_{id}.xyz'),
        png = join(config['path'], 'output', 'mol', 'molid_{id}_2D.png')
    run:
        inchi = read_string(input[0])
        inchi = inchi2geom(inchi, output[0], output[1], output[2], ffield='gaff')

rule calculatepKa:
    input:
        rules.inchi2geom.output.mol3D
    output:
        join(config['path'], 'output', 'pKa', 'molid_{id}_3D.pka')
    shell:
        'cxcalc pka -i -40 -x 40 -d large {input} > {output}'

# rule generateAdducts:
#     input:
#         geometry = rules.InChI2xyz.output.mol_3d,
#         pKa_values = rules.calculatepKa.output[:]
#     output:
#         # sodiated_xyz
#         # protonated_xyz
#         # deprotonated_xyz
#         # sodiated_mol2
#         # protonated_mol2
#         # deprotonated_mol2
#     run:
#         # resources/adduct_creation/create_adducts.py

#         # Final output once we have CCS
#         # InChI, Processed InChI, Formula, Mass, CCS
