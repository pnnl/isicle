from os.path import *
import pybel
import pandas as pd

# snakemake configuration
include: 'mobility.snakefile'

# localrules: all, inchi2smiles, canonicalize, desalt, neutralize, tautomerize,
#             calculateFormula, calculateMass, generateGeometry, generateAdduct,
#             prepare, restore, tleapConfig, sanderEMconfig, sander0config, sanderConfig,
#             selectFrames, extractFrames, convert, calculate_rmsd, downselect, copyOver,
#             createDFTConfig, parseDFT, parseMobcal, boltzmannAverage, calibrate

SMI, = glob_wildcards(abspath(join('input', '{id}.smi')))
INCHI, = glob_wildcards(abspath(join('input', '{id}.inchi')))
IDS = SMI + INCHI

mass = []
for i in IDS:
    mass.append(pybel.readfile(splitext(i)[-1][1:], i).next().molwt)

df = pd.DataFrame({'id': IDS, 'mass': mass})
df.sort_values(by='mass', inplace=True)

IDS = df['id'].values

if 'stop' in config:
    IDS = IDS[config['start']:config['stop']]


rule all:
    input:
        expand(abspath(join('output', 'mobility', 'mobcal', 'calibrated_ccs', '{id}_{adduct}.tsv')),
               id=IDS, adduct=config['adducts'])
