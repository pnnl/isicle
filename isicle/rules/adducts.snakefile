from os.path import *

# snakemake configuration
ruleorder: canonicalize > inchi2smiles


rule inchi2smiles:
    input:
        abspath(join('input', '{id}.inchi'))
    output:
        abspath(join('output', 'adducts', 'canonicalized', '{id}.smi'))
    version:
        "obabel 2> /dev/null | grep 'Open Babel' | awk '{print $3}'"
    log:
        abspath(join('output', 'adducts', 'canonicalized', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'canonicalized', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
        'echo `cat {input}` | obabel -iinchi -ocan > {output} 2> {log}'


rule canonicalize:
    input:
        abspath(join('input', '{id}.smi'))
    output:
        abspath(join('output', 'adducts', 'canonicalized', '{id}.smi'))
    version:
        "obabel 2> /dev/null | grep 'Open Babel' | awk '{print $3}'"
    log:
        abspath(join('output', 'adducts', 'canonicalized', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'canonicalized', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
        'echo `cat {input}` | obabel -ismi -ocan > {output} 2> {log}'


rule desalt:
    input:
        abspath(join('output', 'adducts', 'canonicalized', '{id}.smi'))
    output:
        abspath(join('output', 'adducts', 'desalted', '{id}.smi'))
    version:
        "obabel 2> /dev/null | grep 'Open Babel' | awk '{print $3}'"
    log:
        abspath(join('output', 'adducts', 'desalted', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'desalted', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
        'echo `cat {input}` | obabel -ican -r -ocan > {output} 2> {log}'


rule neutralize:
    input:
        rules.desalt.output
    output:
        abspath(join('output', 'adducts', 'neutralized', '{id}.smi'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'adducts', 'neutralized', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'neutralized', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.process_smiles {input} {output} --neutralize &> {log}'


rule tautomerize:
    input:
        rules.neutralize.output
    output:
        abspath(join('output', 'adducts', 'tautomer', '{id}.smi'))
    log:
        abspath(join('output', 'adducts', 'tautomer', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'tautomer', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
        'obtautomer -c {input} > {output} 2> {log}'


rule calculateFormula:
    input:
        rules.tautomerize.output
    output:
        abspath(join('output', 'adducts', 'formula', '{id}.formula'))
    version:
        "obabel 2> /dev/null | grep 'Open Babel' | awk '{print $3}'"
    log:
        abspath(join('output', 'adducts', 'formula', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'formula', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
       'obabel {input} -o smiles --append "formula" | cut -f2 > {output} 2> {log}'


rule calculateMass:
    input:
        rules.calculateFormula.output
    output:
        abspath(join('output', 'adducts', 'mass', '{id}.mass'))
    version:
        "python -m isicle.scripts.molmass --version | awk '{print $2}'"
    log:
        abspath(join('output', 'adducts', 'mass', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'mass', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.molmass `cat {input}` > {output} 2> {log}'


rule generateGeometry:
    input:
        rules.tautomerize.output
    output:
        mol = abspath(join('output', 'adducts', 'geometry_parent', '{id}.mol')),
        mol2 = abspath(join('output', 'adducts', 'geometry_parent', '{id}.mol2')),
        xyz = abspath(join('output', 'adducts', 'geometry_parent', '{id}.xyz')),
        png = abspath(join('output', 'adducts', 'geometry_parent', 'images', '{id}.png'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'adducts', 'geometry_parent', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'geometry_parent', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.generate_geometry {input} {output.mol} {output.mol2} {output.xyz} {output.png} \
         --forcefield {config[forcefield][type]} --steps {config[forcefield][steps]} &> {log}'


checkpoint generateAdducts:
    input:
        mol2 = rules.generateGeometry.output.mol2
    output:
        directory(join('output', 'adducts', 'geometry_{adduct}','{id}'))
    shell:
        'python -m isicle.scripts.generate_adduct {input.mol2} [{wildcards.adduct}] \
        --forcefield {config[forcefield][type]} --steps {config[forcefield][steps]}' # &> {log}'


rule touchAdducts:
    input:
        xyz = abspath(join('output', 'adducts', 'geometry_{adduct}','{id}', '{addID}.xyz')),
        mol2 = abspath(join('output', 'adducts', 'geometry_{adduct}','{id}', '{addID}.mol2')),
        pdb = abspath(join('output', 'adducts', 'geometry_{adduct}','{id}', '{addID}.pdb')),
        charge = abspath(join('output', 'adducts', 'geometry_{adduct}','{id}', '{addID}.charge'))
    output:
        xyz = abspath(join('output', 'touch', 'geometry_{adduct}', '{id}_{adduct}_{addID}.xyz')),
        mol2 = abspath(join('output', 'touch', 'geometry_{adduct}', '{id}_{adduct}_{addID}.mol2')),
        pdb = abspath(join('output', 'touch', 'geometry_{adduct}', '{id}_{adduct}_{addID}.pdb')),
        charge = abspath(join('output', 'touch', 'geometry_{adduct}', '{id}_{adduct}_{addID}.charge'))
    shell:
        """
        cp -f {input.xyz} {output.xyz}
        cp -f {input.mol2} {output.mol2}
        obabel -ipdb {input.pdb} -opdb -O{output.pdb} -xn
        cp -f {input.charge} {output.charge}
        """