from os.path import *


rule desaltInChI:
    input:
        join(config['path'], 'input', '{id}.inchi')
    output:
        join(config['path'], 'output', 'adducts', 'desalted', '{id}.inchi')
    log:
        join(config['path'], 'output', 'adducts', 'desalted', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'adducts', 'desalted', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python isicle/process_inchi.py {input} {output} --desalt &> {log}'


rule neutralizeInChI:
    input:
        rules.desaltInChI.output
    output:
        join(config['path'], 'output', 'adducts', 'neutralized', '{id}.inchi')
    log:
        join(config['path'], 'output', 'adducts', 'neutralized', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'adducts', 'neutralized', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python isicle/process_inchi.py {input} {output} --neutralize &> {log}'


rule tautomerizeInChI:
    input:
        rules.neutralizeInChI.output
    output:
        join(config['path'], 'output', 'adducts', 'tautomer', '{id}.inchi')
    log:
        join(config['path'], 'output', 'adducts', 'tautomer', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'adducts', 'tautomer', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python isicle/process_inchi.py {input} {output} --tautomerize &> {log}'


rule calculateFormula:
    input:
        rules.tautomerizeInChI.output
    output:
        join(config['path'], 'output', 'adducts', 'formula', '{id}.formula')
    log:
        join(config['path'], 'output', 'adducts', 'formula', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'adducts', 'formula', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python isicle/process_inchi.py {input} {output} --formula &> {log}'


rule calculateMass:
    input:
        rules.calculateFormula.output
    output:
        join(config['path'], 'output', 'adducts', 'mass', '{id}.mass')
    benchmark:
        join(config['path'], 'output', 'adducts', 'mass', 'benchmarks', '{id}.benchmark')
    group:
        'adducts'
    shell:
        'python isicle/resources/molmass.py `cat {input}` &> {output}'


rule generateGeometry:
    input:
        rules.tautomerizeInChI.output
    output:
        mol = join(config['path'], 'output', 'adducts', 'geometry_parent', '{id}.mol'),
        mol2 = join(config['path'], 'output', 'adducts', 'geometry_parent', '{id}.mol2'),
        xyz = join(config['path'], 'output', 'adducts', 'geometry_parent', '{id}.xyz'),
        png = join(config['path'], 'output', 'adducts', 'geometry_parent', 'images', '{id}.png')
    log:
        join(config['path'], 'output', 'adducts', 'geometry_parent', 'logs', '{id}.log')
    benchmark:
        join(config['path'], 'output', 'adducts', 'geometry_parent', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python isicle/generate_geometry.py {input} {output.mol} {output.mol2} {output.xyz} {output.png} \
         --forcefield {config[forcefield][type]} --steps {config[forcefield][steps]} &> {log}'


rule calculatepKa:
    input:
        rules.generateGeometry.output.mol
    output:
        join(config['path'], 'output', 'adducts', 'pKa', '{id}.pka')
    benchmark:
        join(config['path'], 'output', 'adducts', 'pKa', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'cxcalc pka -i -40 -x 40 -d large {input} &> {output}'


rule generateAdduct:
    input:
        molfile = rules.generateGeometry.output.mol,
        pkafile = rules.calculatepKa.output
    output:
        xyz = join(config['path'], 'output', 'adducts', 'geometry_{adduct}', '{id}_{adduct}.xyz'),
        mol2 = join(config['path'], 'output', 'adducts', 'geometry_{adduct}', '{id}_{adduct}.mol2')
    log:
        join(config['path'], 'output', 'adducts', 'geometry_{adduct}', 'logs', '{id}_{adduct}.log')
    benchmark:
        join(config['path'], 'output', 'adducts', 'geometry_{adduct}', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python isicle/generate_adduct.py {input.molfile} {input.pkafile} [{wildcards.adduct}] {output.mol2} {output.xyz} \
         --forcefield {config[forcefield][type]} --steps {config[forcefield][steps]} &> {log}'
