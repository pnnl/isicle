from os.path import *
from pkg_resources import resource_filename

# snakemake configuration
include: 'molecular_dynamics.snakefile', 'adducts.snakefile'


rule copyOver:
    input:
        abspath(join('output', 'md', 'downselected', '{id}_{adduct}_{addID}_{cycle}_{selected}.xyz'))
    output:
        abspath(join('output', 'dft', '{id}_{adduct}_{addID}', 'cycle_{cycle}_{selected}', '{id}_{adduct}_{addID}_{cycle}_{selected}.xyz'))
    log:
        abspath(join('output', 'dft', 'logs', '{id}_{adduct}_{addID}_{cycle}_{selected}.copy.log'))
    benchmark:
        abspath(join('output', 'dft', 'benchmarks', '{id}_{adduct}_{addID}_{cycle}_{selected}.copy.benchmark'))
    # group:
    #     'dft'
    shell:
        'cp {input} {output} &> {log}'


# create .nw files based on template (resources/nwchem/template.nw)
rule createDFTConfig:
    input:
        xyz = rules.copyOver.output,
        charge = rules.touchAdducts.output.charge
    output:
        abspath(join('output', 'dft', '{id}_{adduct}_{addID}', 'cycle_{cycle}_{selected}', '{id}_{adduct}_{addID}_{cycle}_{selected}.nw'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'dft', 'logs', '{id}_{adduct}_{addID}_{cycle}_{selected}.create.log'))
    benchmark:
        abspath(join('output', 'dft', 'benchmarks', '{id}_{adduct}_{addID}_{cycle}_{selected}.create.benchmark'))
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
        abspath(join('output', 'dft', '{id}_{adduct}_{addID}', 'cycle_{cycle}_{selected}', '{id}_{adduct}_{addID}_{cycle}_{selected}.out'))
    version:
        "nwchem /dev/null | grep '(NWChem)' | awk '{print $6}'"
    log:
        abspath(join('output', 'dft', 'logs', '{id}_{adduct}_{addID}_{cycle}_{selected}.nwchem.log'))
    benchmark:
        abspath(join('output', 'dft', 'benchmarks', '{id}_{adduct}_{addID}_{cycle}_{selected}.nwchem.benchmark'))
    # group:
    #     'dft'
    threads: 8
    shell:
        'nwchem {input} > {output} 2> {log}'


# parse nwchem outputs (geometry files)
rule parseDFT:
    input:
        rules.dft.output
    output:
        mfj = abspath(join('output', 'mobility', 'mobcal', 'runs', '{id}_{adduct}_{addID}_{cycle}_{selected}.mfj')),
        energy = abspath(join('output', 'mobility', 'mobcal', 'runs', '{id}_{adduct}_{addID}_{cycle}_{selected}.energy'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'dft', 'logs', '{id}_{adduct}_{addID}_{cycle}_{selected}.parse.log'))
    benchmark:
        abspath(join('output', 'dft', 'benchmarks', '{id}_{adduct}_{addID}_{cycle}_{selected}.parse.benchmark'))
    # group:
    #     'dft'
    shell:
        'python -m isicle.scripts.parse_nwchem {input} {output.mfj} --dft &> {log}'
