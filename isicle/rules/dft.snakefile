from os.path import *

# snakemake configuration
include: 'molecular_dynamics.snakefile'


rule copyOver:
    input:
        join(config['path'], 'output', 'md', 'downselected', '{id}_{adduct}_{cycle}_{selected}.xyz')
    output:
        join(config['path'], 'output', 'dft', '{id}_{adduct}', 'cycle_{cycle}_{selected}', '{id}_{adduct}_{cycle}_{selected}.xyz')
    # group:
    #     'dft'
    shell:
        'cp {input} {output}'


# create .nw files based on template (resources/nwchem/template.nw)
rule createNW:
    input:
        rules.copyOver.output
    output:
        join(config['path'], 'output', 'dft', '{id}_{adduct}', 'cycle_{cycle}_{selected}', '{id}_{adduct}_{cycle}_{selected}.nw')
    log:
        join(config['path'], 'output', 'dft', 'logs', '{id}_{adduct}_{cycle}_{selected}.create.log')
    benchmark:
        join(config['path'], 'output', 'dft', 'benchmarks', '{id}_{adduct}_{cycle}_{selected}.create.benchmark')
    # group:
    #     'dft'
    shell:
        'python isicle/generateNW.py {input} --template {config[nwchem][dft_template]} > {log}'


# run NWChem
rule NWChem:
    input:
        rules.createNW.output
    output:
        join(config['path'], 'output', 'dft', '{id}_{adduct}', 'cycle_{cycle}_{selected}', '{id}_{adduct}_{cycle}_{selected}.out')
    benchmark:
        join(config['path'], 'output', 'dft', 'benchmarks', '{id}_{adduct}_{cycle}_{selected}.dft.benchmark')
    # group:
    #     'dft'
    shell:
        'srun --mpi=pmi2 nwchem {input} > {output}'


# parse nwchem outputs (geometry files)
rule parseNWChem:
    input:
        rules.NWChem.output
    output:
        geom1 = join(config['path'], 'output', 'mobility', 'mobcal', 'runs', '{id}_{adduct}_{cycle}_{selected}_charge.mfj'),
        geom2 = join(config['path'], 'output', 'mobility', 'mobcal', 'runs', '{id}_{adduct}_{cycle}_{selected}_geom+charge.mfj'),
        charge1 = join(config['path'], 'output', 'mobility', 'mobcal', 'runs', '{id}_{adduct}_{cycle}_{selected}_charge.energy'),
        charge2 = join(config['path'], 'output', 'mobility', 'mobcal', 'runs', '{id}_{adduct}_{cycle}_{selected}_geom+charge.energy')
    log:
        join(config['path'], 'output', 'dft', 'logs', '{id}_{adduct}_{cycle}_{selected}.parse.log')
    benchmark:
        join(config['path'], 'output', 'dft', 'benchmarks', '{id}_{adduct}_{cycle}_{selected}.parse.benchmark')
    # group:
    #     'dft'
    run:
        outdir = dirname(output.geom2)
        shell('python isicle/parse_nwchem.py {input} %s --mode dft > {log}' % outdir)
