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
        inchi = expand(rules.desalt.input, id=IDS),
        p_inchi = expand(rules.tautomerize.output, id=IDS)
    output:
        join(config['path'], 'output', 'impact_results.tsv')
    run:
        dfs = []
        for f in input.ccs:
            # read ccs
            ccs = read_impact(f)

            # additional info
            ID, adduct = splitext(basename(f))[0].rsplit('_', 1)
            mass = read_mass([x for x in input.mass if ID in x][0])
            formula = read_string([x for x in input.formula if ID in x][0])
            inchi = read_string([x for x in input.inchi if ID in x][0])
            p_inchi = read_string([x for x in input.p_inchi if ID in x][0])

            df = pd.DataFrame(data=[[ID, adduct, inchi, formula, mass, p_inchi, ccs]],
                              columns=['ID', 'Adduct', 'Parent InChI', 'Parent Formula',
                                       'Parent Mass', 'Processed InChI', 'CCS_He'])

            # reorder/rename
            df['CCS_N2'] = df['CCS_He'] + config['ccs']['alpha'] * df['Parent Mass'] ** config['ccs']['beta']
            dfs.append(df)

        master = pd.concat(dfs)
        master.to_csv(output[0], sep='\t', index=False)
        # newf = open(output[0], 'w')
        # newf.write(master.to_string(col_space=5, justify='left'))
        # newf.close()
