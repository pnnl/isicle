from os.path import *
from resources.mobcal.finish import boltzmann, parse_mobcal
import pandas as pd
from resources.utils import *

# snakemake configuration
configfile: 'config.yaml'
include: 'dft.snakefile'


IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

# run mobcal on geom+charge nwchem output
rule mobcal:
    input:
        rules.parseNWChem.output.geom2
    output:
        join(config['path'], 'output', 'mobcal', '{id}_{adduct}_{cycle}_{selected}_geom+charge.out')
    group:
        'mobility'
    shell:
        '{config[mobcal][runscript]} {config[mobcal][params]} {config[mobcal][atomtypes]} {input}'

# parse mobcal output
rule parseMobcal:
    input:
        geom = expand(rules.mobcal.output, id=IDS, adduct=config['adducts'], cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2']),
        energy = expand(rules.parseNWChem.output.charge2, id=IDS, adduct=config['adducts'], cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2'])
    output:
        join(config['path'], 'output', 'ccs_all_conformers.tsv')
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
            else:
                print(ccsfile)

        df = pd.DataFrame(res, columns=['File', 'Mobility', 'Mean CCS', 'Stdev CCS', 'DFT Energy'])
        df.to_csv(output[0], sep='\t', index=False)

# boltzmann averaging
rule boltzmann_average:
    input:
        rules.parseMobcal.output
    output:
        join(config['path'], 'output', 'ccs_result.tsv')
    group:
        'mobility'
    run:
        boltzmann(input[0], output[0])
