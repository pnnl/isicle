from os.path import *
import pandas as pd
from resources.utils import *

# snakemake configuration
configfile: 'config.yaml'
localrules: impact, postprocess
include: 'adducts.snakefile'

# infer wildcards from inputs, select adducts
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))
OS = getOS()

rule impact:
    input:
        rules.generateAdducts.output.xyz
    output:
        join(config['path'], 'output', '5_impact', '{id}_{adduct}.txt')
    benchmark:
        join(config['path'], 'output', '5_impact', 'benchmarks', '{id}_{adduct}.benchmark')
    shell:
        # run impact on adducts
        'resources/IMPACT/{OS}/impact {input} -o {output} -H \
        -shotsPerRot {config[impact][shotsPerRot]} \
        -convergence {config[impact][convergence]} \
        -nRuns {config[impact][nRuns]} -nocite'

rule postprocess:
    input:
        ccs = expand(rules.impact.output, id=IDS, adduct=config['adducts']),
        mass = expand(rules.calculateMass.output, id=IDS),
        formula = expand(rules.calculateFormula.output, id=IDS),
        inchi = expand(rules.desalt.input, id=IDS)
    output:
        join(config['path'], 'output', 'impact_results.tsv')
    run:
        dfs = []
        for ccs in input.ccs:
            # read ccs
            df = read_impact(ccs)

            # ID, adduct
            ID, adduct = splitext(basename(ccs))[0].rsplit('_', 1)

            # extract
            mass = [x for x in input.mass if ID in x][0]
            formula = [x for x in input.formula if ID in x][0]
            inchi = [x for x in input.inchi if ID in x][0]

            # additional info
            df['ID'] = ID
            df['Adduct'] = adduct
            df['Parent Mass'] = read_mass(mass)
            df['Parent Formula'] = read_string(formula)
            df['Parent InChI'] = read_string(inchi)

            # reorder/rename
            df = df[['ID', 'Adduct', 'Parent InChI', 'Parent Formula', 'Parent Mass', 'CCS_TJM']]
            df.columns = ['ID', 'Adduct', 'Parent InChI', 'Parent Formula', 'Parent Mass', 'CCS_He']
            df['CCS_N2'] = df['CCS_He'] + config['ccs']['alpha'] * df['Parent Mass'] ** config['ccs']['beta']
            dfs.append(df)

        master = pd.concat(dfs)
        master.to_csv(output[0], sep='\t', index=False)
        # newf = open(output[0], 'w')
        # newf.write(master.to_string(col_space=5, justify='left'))
        # newf.close()
