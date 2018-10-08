from os.path import *
import pandas as pd
from resources.utils import *

# snakemake configuration
configfile: 'config.yaml'
include: 'adducts.snakefile'

OS = getOS()


rule impact:
    input:
        rules.generateAdducts.output.xyz
    output:
        join(config['path'], 'output', 'mobility', 'impact', 'runs', '{id}_{adduct}.txt')
    benchmark:
        join(config['path'], 'output', 'mobility', 'impact', 'runs', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'mobility_alt'
    shell:
        # run impact on adducts
        'IMPACT_RANDSEED={config[impact][seed]} resources/IMPACT/{OS}/impact {input} -o {output} -H \
        -shotsPerRot {config[impact][shotsPerRot]} \
        -convergence {config[impact][convergence]} \
        -nRuns {config[impact][nRuns]} -nocite'

rule postprocess:
    input:
        ccs = rules.impact.output,
        mass = rules.calculateMass.output
    output:
        join(config['path'], 'output', 'mobility', 'impact', 'ccs', '{id}_{adduct}.N2.ccs'),
        join(config['path'], 'output', 'mobility', 'impact', 'ccs', '{id}_{adduct}.He.ccs')
    benchmark:
        join(config['path'], 'output', 'mobility', 'impact', 'ccs', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'mobility_alt'
    run:
        # read inputs
        ccs_He = read_impact(input[0])
        m = read_mass(input[1])

        ccs_N2 = ccs_He + config['ccs']['alpha'] * m ** config['ccs']['beta']

        write_string(str(ccs_N2), output[0])
        write_string(str(ccs_He), output[1])

# # for report
# ID, adduct = splitext(basename(f))[0].rsplit('_', 1)
# mass = read_mass([x for x in input.mass if ID in x][0])
# formula = read_string([x for x in input.formula if ID in x][0])
# inchi = read_string([x for x in input.inchi if ID in x][0])
# p_inchi = read_string([x for x in input.p_inchi if ID in x][0])

# df = pd.DataFrame(data=[[ID, adduct, inchi, formula, mass, p_inchi, ccs]],
#                   columns=['ID', 'Adduct', 'Parent InChI', 'Parent Formula',
#                            'Parent Mass', 'Processed InChI', 'CCS_He'])
