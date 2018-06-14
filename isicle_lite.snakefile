from os.path import *
from glob import glob

# snakemake configuration
configfile: 'config.yaml'
localrules: all, impact, process
include: 'adducts.snakefile'

# infer wildcards from inputs, select adducts
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))
ADDUCTS = ['+Na', '+H', '-H']

# a pseudo-rule that collects the target files
rule all:
    input:
        join(config['path'], 'impact.txt')

rule impact:
    input:
        expand(join(config['path'], 'output', 'adduct_structures', 'xyz', '{id}_{adduct}.xyz'), id=IDS, adduct=ADDUCTS)
    output:
        join(config['path'], 'impact.txt')
    shell:
        # run impact exe on adduct folder
        'resources/IMPACT/impact {input} -o impact.txt -H -shotsPerRot 64 -convergence .001 -nRuns 64'

rule process:
    input:
        rules.impact.output
    output:
        join(config['path'], 'impact_processed.txt')
    run:
        with open(input[0]) as f:
            flines = f.readlines()[8:]

        # Process and save in processed.txt file
        with open(output[0], 'w') as f:
            for line in flines:
                tokens = [token for token in line.split(' ') if len(token) > 0]
                molid = line.split('\\')[-1].split('+')[0]
                ccs = tokens[-2]
                f.write('\t'.join([molid, ccs]))
                f.write('\n')
