from os.path import *
from resources.mobcal.finish import boltzmann, parse_mobcal
import pandas as pd
from resources.utils import *

# snakemake configuration
configfile: 'isicle/config.yaml'
include: 'dft.snakefile'


# run mobcal on geom+charge nwchem output
rule mobcal:
    input:
        rules.parseNWChem.output.geom2
    output:
        join(config['path'], 'output', 'mobcal', '{id}_{adduct}_{cycle}_{selected}_geom+charge.out')
    benchmark:
        join(config['path'], 'output', 'mobcal', 'benchmarks', '{id}_{adduct}_{cycle}_{selected}.benchmark')
    group:
        'mobility'
    shell:
        '{config[mobcal][runscript]} {config[mobcal][params]} {config[mobcal][atomtypes]} {input}'

# parse mobcal output
rule parseMobcal:
    input:
        geom = expand(join(config['path'], 'output', 'mobcal', '{{id}}_{{adduct}}_{cycle}_{selected}_geom+charge.out'),
                      cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2']),
        energy = expand(join(config['path'], 'output', 'mobcal', '{{id}}_{{adduct}}_{cycle}_{selected}_geom+charge.energy'),
                        cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2'])
    output:
        join(config['path'], 'output', 'conformer_ccs', '{id}_{adduct}.tsv')
    benchmark:
        join(config['path'], 'output', 'conformer_ccs', 'benchmarks', '{id}_{adduct}.benchmark')
    group:
        'mobility'
    run:
        res = []
        for ccsfile, efile in zip(input['geom'], input['energy']):
            tmp = parse_mobcal(ccsfile)
            if tmp is not None:
                # get energy
                with open(efile, 'r') as f:
                    e = float(f.readlines()[0])
                tmp.append(e)
                res.append(tmp)

        df = pd.DataFrame(res, columns=['Mobility', 'Mean CCS', 'Stdev CCS', 'DFT Energy'])
        df.to_csv(output[0], sep='\t', index=False)

# boltzmann averaging
rule boltzmannAverage:
    input:
        rules.parseMobcal.output
    output:
        join(config['path'], 'output', 'boltzmann_ccs', '{id}_{adduct}.tsv')
    benchmark:
        join(config['path'], 'output', 'boltzmann_ccs', 'benchmarks', '{id}_{adduct}.benchmark')
    group:
        'mobility'
    run:
        boltzmann(input[0], output[0])
