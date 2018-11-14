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
    version:
        "cxcalc --help | grep 'version ' | awk '{print $2}'"
    log:
        abspath(join('output', 'adducts', 'tautomer', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'tautomer', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
        'cxcalc majortautomer -f smiles `cat {input}` > {output} 2> {log}'


rule calculateFormula:
    input:
        rules.tautomerize.output
    output:
        abspath(join('output', 'adducts', 'formula', '{id}.formula'))
    version:
        "cxcalc --help | grep 'version ' | awk '{print $2}'"
    log:
        abspath(join('output', 'adducts', 'formula', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'formula', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
        "cxcalc formula `cat {input}` | tail -n1 | awk '{{print $2}}' > {output} 2> {log}"


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


rule calculatepKa:
    input:
        rules.generateGeometry.output.mol
    output:
        abspath(join('output', 'adducts', 'pKa', '{id}.pka'))
    version:
        "cxcalc --help | grep 'version ' | awk '{print $2}'"
    log:
        abspath(join('output', 'adducts', 'pKa', 'logs', '{id}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'pKa', 'benchmarks', '{id}.benchmark'))
    # group:
    #     'adducts'
    shell:
        'cxcalc pka -i -40 -x 40 -d large {input} > {output} 2> {log}'


rule generateAdduct:
    input:
        molfile = rules.generateGeometry.output.mol,
        pkafile = rules.calculatepKa.output
    output:
        xyz = abspath(join('output', 'adducts', 'geometry_{adduct}', '{id}_{adduct}.xyz')),
        mol2 = abspath(join('output', 'adducts', 'geometry_{adduct}', '{id}_{adduct}.mol2'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'adducts', 'geometry_{adduct}', 'logs', '{id}_{adduct}.log'))
    benchmark:
        abspath(join('output', 'adducts', 'geometry_{adduct}', 'benchmarks', '{id}_{adduct}.benchmark'))
    # group:
    #     'adducts'
    shell:
        'python -m isicle.scripts.generate_adduct {input.molfile} {input.pkafile} [{wildcards.adduct}] {output.mol2} {output.xyz} \
         --forcefield {config[forcefield][type]} --steps {config[forcefield][steps]} &> {log}'
