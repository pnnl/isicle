from os.path import *
from resources.utils import *
import sys
import logging

# snakemake configuration
configfile: 'isicle/config.yaml'


rule desaltInChI:
    input:
        join(config['path'], 'input', '{id}.inchi')
    output:
        join(config['path'], 'output', '{id}', '0_desalted', '{id}.inchi')
    log:
        join(config['path'], 'output', '{id}', '0_desalted', '{id}.log')
    benchmark:
        join(config['path'], 'output', '{id}', '0_desalted', '{id}.benchmark')
    group:
        'adducts'
    run:
        inchi = read_string(input[0])
        inchi = desalt(inchi, log[0])
        if inchi is not None:
            write_string(inchi, output[0])
        else:
            sys.exit(1)

rule neutralizeInChI:
    input:
        rules.desaltInChI.output
    output:
        join(config['path'], 'output', '{id}', '1_neutralized', '{id}.inchi')
    benchmark:
        join(config['path'], 'output', '{id}', '1_neutralized', '{id}.benchmark')
    group:
        'adducts'
    run:
        inchi = read_string(input[0])
        inchi = neutralize(inchi)
        if inchi is not None:
            write_string(inchi, output[0])
        else:
            sys.exit(1)

rule tautomerizeInChI:
    input:
        rules.neutralizeInChI.output
    output:
        join(config['path'], 'output', '{id}', '2a_tautomer', '{id}.inchi')
    log:
        join(config['path'], 'output', '{id}', '2a_tautomer', '{id}.log')
    benchmark:
        join(config['path'], 'output', '{id}', '2a_tautomer', '{id}.benchmark')
    group:
        'adducts'
    run:
        inchi = read_string(input[0])
        inchi = tautomerize(inchi, log=log[0])
        if inchi is not None:
            write_string(inchi, output[0])
        else:
            sys.exit(1)

rule calculateFormula:
    input:
        rules.tautomerizeInChI.output
    output:
        join(config['path'], 'output', '{id}', 'parent', '2b_formula', '{id}.formula')
    log:
        join(config['path'], 'output', '{id}', 'parent', '2b_formula', '{id}.log')
    benchmark:
        join(config['path'], 'output', '{id}', 'parent', '2b_formula', '{id}.benchmark')
    group:
        'adducts'
    run:
        inchi = read_string(input[0])
        formula = inchi2formula(inchi, log=log[0])
        write_string(formula, output[0])

rule calculateMass:
    input:
        rules.calculateFormula.output
    output:
        join(config['path'], 'output', '{id}', 'parent', '2b_mass', '{id}.mass')
    benchmark:
        join(config['path'], 'output', '{id}', 'parent', '2b_mass', '{id}.benchmark')
    group:
        'adducts'
    shell:
        'python isicle/resources/molmass.py `cat {input}` > {output}'

rule generateGeometry:
    input:
        rules.tautomerizeInChI.output
    output:
        mol = join(config['path'], 'output', '{id}', 'parent', '3a_geometry', '{id}.mol'),
        png = join(config['path'], 'output', '{id}', 'parent', '3b_image', '{id}.png')
    benchmark:
        join(config['path'], 'output', '{id}', 'parent', '3a_geometry', '{id}.benchmark')
    group:
        'adducts'
    run:
        inchi = read_string(input[0])
        mol = inchi2geom(inchi, forcefield=config['forcefield']['type'],
                         steps=config['forcefield']['steps'])

        mol.draw(show=False, filename=output.png, usecoords=False, update=False)
        mol.write('mol', output.mol, overwrite=True)

rule calculatepKa:
    input:
        rules.generateGeometry.output.mol
    output:
        join(config['path'], 'output', '{id}', 'parent', '3c_pKa', '{id}.pka')
    benchmark:
        join(config['path'], 'output', '{id}', 'parent', '3c_pKa', '{id}.benchmark')
    group:
        'adducts'
    shell:
        'cxcalc pka -i -40 -x 40 -d large {input} > {output}'

rule generateAdducts:
    input:
        molfile = rules.generateGeometry.output.mol,
        pkafile = rules.calculatepKa.output
    output:
        xyz = join(config['path'], 'output', '{id}', 'adduct_{adduct}', '0_geometry', '{id}_{adduct}.xyz'),
        mol2 = join(config['path'], 'output', '{id}', 'adduct_{adduct}', '0_geometry', '{id}_{adduct}.mol2')
    log:
        join(config['path'], 'output', '{id}', 'adduct_{adduct}', '0_geometry', '{id}_{adduct}.log')
    benchmark:
        join(config['path'], 'output', '{id}', 'adduct_{adduct}', '0_geometry', '{id}_{adduct}.benchmark')
    group:
        'adducts'
    run:
        # log
        logging.basicConfig(filename=log[0], level=logging.DEBUG)

        # read inputs
        mol = next(pybel.readfile("mol", input[0]))
        pka = read_pka(input[1])

        # generate adduct
        a = wildcards.adduct
        if '+' in a:
            adduct = create_adduct(mol, a, pka['b1'],
                                   forcefield=config['forcefield']['type'],
                                   steps=config['forcefield']['steps'])
        elif '-' in a:
            adduct = create_adduct(mol, a, pka['a1'],
                                   forcefield=config['forcefield']['type'],
                                   steps=config['forcefield']['steps'])
        else:
            adduct = mol

        adduct.write('xyz', output.xyz, overwrite=True)
        adduct.write('mol2', output.mol2, overwrite=True)
