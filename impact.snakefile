from os.path import *
import pandas as pd
from resources.utils import *

# snakemake configuration
configfile: 'config.yaml'
localrules: impact, postprocess
include: 'adducts.snakefile'

# infer wildcards from inputs, select adducts
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

rule impact:
    input:
        rules.generateAdducts.output.xyz
    output:
        join(config['path'], 'output', '5_impact', '{id}_{adduct}.txt')
    shell:
        # run impact on adducts
        'resources/IMPACT/impact {input} -o {output} -H -shotsPerRot 64 -convergence .001 -nRuns 64 -nocite'

rule postprocess:
    input:
        ccs = expand(rules.impact.output, id=IDS, adduct=config['adducts']),
        mass = expand(rules.calculateMass.output, id=IDS),
        formula = expand(rules.calculateFormula.output, id=IDS),
        inchi = expand(rules.tautomerize.output, id=IDS)
    output:
        join(config['path'], 'output', 'impact_results.tsv')
    run:
        dfs = []
        for ccs, mass, formula, inchi, ID in zip(input.ccs, input.mass,
                                                 input.formula, input.inchi, IDS):
            # read ccs
            df = read_impact(ccs)

            # additional info
            df['ID'] = ID
            df['Adduct'] = splitext(basename(ccs))[0].split('_')[-1]
            df['Parent Mass'] = read_mass(mass)
            df['Parent Formula'] = read_string(formula)
            df['Parent InChI'] = read_string(inchi)

            # reorder
            df = df[['ID', 'Adduct', 'Parent InChI', 'Parent Formula', 'Parent Mass', 'CCS_PA', 'SEM_rel', 'CCS_TJM']]
            dfs.append(df)

        master = pd.concat(dfs)
        master.to_csv(output[0], sep='\t', index=False)
