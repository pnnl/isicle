from os.path import *
import pandas as pd

# snakemake configuration
configfile: 'config.yaml'
localrules: all, impact, process
include: 'adducts.snakefile'

rule impact:
    input:
        rules.generateAdducts.output.xyz
    output:
        join(config['path'], 'output', '5_impact', '{id}_{adduct}.txt')
    shell:
        # run impact on adducts
        'resources/IMPACT/impact {input} -o {output} -H -shotsPerRot 64 -convergence .001 -nRuns 64 -nocite'

# rule process:
#     input:
#         rules.impact.output
#     output:
#         join(config['path'], 'impact.txt')
#     run:
#         df = pd.read_csv(input[0], sep='\t')
