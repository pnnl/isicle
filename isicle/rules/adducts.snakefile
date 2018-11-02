from os.path import *

# snakemake configuration
ruleorder: canonicalize > inchi2smiles


rule inchi2smiles:
    input:
        join('input', '{id}.inchi')
    output:
        join('output', 'adducts', 'canonicalized', '{id}.smi')
    version:
        'python -m isicle.scripts.process_smiles --version'
    log:
        join('output', 'adducts', 'canonicalized', 'logs', '{id}.log')
    benchmark:
        join('output', 'adducts', 'canonicalized', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.process_smiles {input} {output} --inchi &> {log}'


rule canonicalize:
    input:
        join('input', '{id}.smi')
    output:
        join('output', 'adducts', 'canonicalized', '{id}.smi')
    version:
        'python -m isicle.scripts.process_smiles --version'
    log:
        join('output', 'adducts', 'canonicalized', 'logs', '{id}.log')
    benchmark:
        join('output', 'adducts', 'canonicalized', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.process_smiles {input} {output} --canonicalize &> {log}'


rule desalt:
    input:
        join('output', 'adducts', 'canonicalized', '{id}.smi')
    output:
        join('output', 'adducts', 'desalted', '{id}.smi')
    version:
        'python -m isicle.scripts.process_smiles --version'
    log:
        join('output', 'adducts', 'desalted', 'logs', '{id}.log')
    benchmark:
        join('output', 'adducts', 'desalted', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.process_smiles {input} {output} --desalt &> {log}'


rule neutralize:
    input:
        rules.desalt.output
    output:
        join('output', 'adducts', 'neutralized', '{id}.smi')
    version:
        'python -m isicle.scripts.process_smiles --version'
    log:
        join('output', 'adducts', 'neutralized', 'logs', '{id}.log')
    benchmark:
        join('output', 'adducts', 'neutralized', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.process_smiles {input} {output} --neutralize &> {log}'


rule tautomerize:
    input:
        rules.neutralize.output
    output:
        join('output', 'adducts', 'tautomer', '{id}.smi')
    version:
        'python -m isicle.scripts.process_smiles --version'
    log:
        join('output', 'adducts', 'tautomer', 'logs', '{id}.log')
    benchmark:
        join('output', 'adducts', 'tautomer', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.process_smiles {input} {output} --tautomerize &> {log}'


rule calculateFormula:
    input:
        rules.tautomerize.output
    output:
        join('output', 'adducts', 'formula', '{id}.formula')
    version:
        'python -m isicle.scripts.process_smiles --version'
    log:
        join('output', 'adducts', 'formula', 'logs', '{id}.log')
    benchmark:
        join('output', 'adducts', 'formula', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.process_smiles {input} {output} --formula &> {log}'


rule calculateMass:
    input:
        rules.calculateFormula.output
    output:
        join('output', 'adducts', 'mass', '{id}.mass')
    version:
        'python -m isicle.scripts.molmass --version'
    log:
        join('output', 'adducts', 'mass', 'logs', '{id}.log')
    benchmark:
        join('output', 'adducts', 'mass', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.molmass `cat {input}` > {output} 2> {log}'


rule generateGeometry:
    input:
        rules.tautomerize.output
    output:
        mol = join('output', 'adducts', 'geometry_parent', '{id}.mol'),
        mol2 = join('output', 'adducts', 'geometry_parent', '{id}.mol2'),
        xyz = join('output', 'adducts', 'geometry_parent', '{id}.xyz'),
        png = join('output', 'adducts', 'geometry_parent', 'images', '{id}.png')
    version:
        'python -m isicle.scripts.generate_geometry --version'
    log:
        join('output', 'adducts', 'geometry_parent', 'logs', '{id}.log')
    benchmark:
        join('output', 'adducts', 'geometry_parent', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.generate_geometry {input} {output.mol} {output.mol2} {output.xyz} {output.png} \
         --forcefield {config[forcefield][type]} --steps {config[forcefield][steps]} &> {log}'


rule calculatepKa:
    input:
        rules.generateGeometry.output.mol
    output:
        join('output', 'adducts', 'pKa', '{id}.pka')
    version:
        "cxcalc --help | grep 'version ' | awk '{print $2}'"
    log:
        join('output', 'adducts', 'pKa', 'logs', '{id}.log')
    benchmark:
        join('output', 'adducts', 'pKa', 'benchmarks', '{id}.benchmark')
    # group:
    #     'adducts'
    shell:
        'cxcalc pka -i -40 -x 40 -d large {input} > {output} 2> {log}'


rule generateAdduct:
    input:
        molfile = rules.generateGeometry.output.mol,
        pkafile = rules.calculatepKa.output
    output:
        xyz = join('output', 'adducts', 'geometry_{adduct}', '{id}_{adduct}.xyz'),
        mol2 = join('output', 'adducts', 'geometry_{adduct}', '{id}_{adduct}.mol2')
    version:
        'python -m isicle.scripts.generate_adduct --version'
    log:
        join('output', 'adducts', 'geometry_{adduct}', 'logs', '{id}_{adduct}.log')
    benchmark:
        join('output', 'adducts', 'geometry_{adduct}', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.generate_adduct {input.molfile} {input.pkafile} [{wildcards.adduct}] {output.mol2} {output.xyz} \
         --forcefield {config[forcefield][type]} --steps {config[forcefield][steps]} &> {log}'
