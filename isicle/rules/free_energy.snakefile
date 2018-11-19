from os.path import *

localrules: all, inchi2smiles, canonicalize, generateGeometry,
            copyOver, createDFTConfig, parseDFT

SMI, = glob_wildcards(abspath(join('input', '{id}.smi')))
INCHI, = glob_wildcards(abspath(join('input', '{id}.inchi')))
IDS = SMI + INCHI

IDS.sort()

if 'stop' in config:
    IDS = IDS[config['start']:config['stop']]


rule all:
    input:
        expand(abspath(join('output', 'energy', '{id}.out')), id=IDS)


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


rule generateGeometry:
    input:
        rules.canonicalize.output
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
    run:
        shell('obabel {input} -O {output.mol2} --gen3d --ff {config[forcefield][type]} --steps {config[forcefield][steps]} \
               --partialcharge eem &> {log}')
        shell('obabel {output.mol2} -O {output.mol} >> {log} 2>&1')
        shell('obabel {output.mol2} -O {output.xyz} >> {log} 2>&1')
        shell('obabel {input} -O {output.png} >> {log} 2>&1')

rule charge:
    input:
        canonicalize.rules.output
    output:
        abspath(join('output', 'adducts', 'geometry_parent', '{id}.charge'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'adducts', 'geometry_parent', 'logs', '{id}.charge.log'))
    benchmark:
        abspath(join('output', 'adducts', 'geometry_parent', 'benchmarks', '{id}.charge.benchmark'))
    shell:
        "cxcalc formalcharge {input} | tail -n1 | awk '{print $2}' > {output} 2> {log}"


rule copyOver:
    input:
        rules.generateGeometry.output.xyz
    output:
        abspath(join('output', 'dft', '{id}', '{id}.xyz'))
    log:
        abspath(join('output', 'dft', 'logs', '{id}.copy.log'))
    benchmark:
        abspath(join('output', 'dft', 'benchmarks', '{id}.copy.benchmark'))
    # group:
    #     'dft'
    shell:
        'cp {input} {output} &> {log}'


# create .nw files based on template (resources/nwchem/template.nw)
rule createDFTConfig:
    input:
        xyz = rules.copyOver.output,
        charge = rules.generateAdduct.output.charge
    output:
        abspath(join('output', 'dft', '{id}', '{id}.nw'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'dft', 'logs', '{id}.create.log'))
    benchmark:
        abspath(join('output', 'dft', 'benchmarks', '{id}.create.benchmark'))
    # group:
    #     'dft'
    shell:
        'python -m isicle.scripts.generateNW {input.xyz} --dft --charge `cat {input.charge}` \
         --template {config[nwchem][dft_template]} &> {log}'


# run NWChem
rule dft:
    input:
        rules.createDFTConfig.output
    output:
        abspath(join('output', 'dft', '{id}', '{id}.out'))
    version:
        "nwchem /dev/null | grep '(NWChem)' | awk '{print $6}'"
    log:
        abspath(join('output', 'dft', 'logs', '{id}.nwchem.log'))
    benchmark:
        abspath(join('output', 'dft', 'benchmarks', '{id}.nwchem.benchmark'))
    # group:
    #     'dft'
    shell:
        'srun --mpi=pmi2 nwchem {input} > {output} 2> {log}'


# parse nwchem outputs (geometry files)
rule parseDFT:
    input:
        rules.dft.output
    output:
        mfj = abspath(join('output', 'dft', '{id}', '{id}.mfj')),
        energy = abspath(join('output', 'energy', '{id}.out'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'dft', 'logs', '{id}.parse.log'))
    benchmark:
        abspath(join('output', 'dft', 'benchmarks', '{id}.parse.benchmark'))
    # group:
    #     'dft'
    shell:
        'python -m isicle.scripts.parse_nwchem {input} {output.mfj} --dft &> {log}'
