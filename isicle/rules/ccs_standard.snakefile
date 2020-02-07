from os.path import *
import pybel
import pandas as pd
import glob

# snakemake configuration
include: 'mobility.snakefile'
localrules: all, inchi2smiles, canonicalize, desalt, neutralize, calculateMass, generateGeometry,
            generateAdduct, prepare, restore, tleapConfig, sanderEMconfig, sander0config, sanderConfig,
            selectFrames, extractFrames, convert, calculate_rmsd, downselect, copyOver,
            createDFTConfig, parseDFT, parseMobcal, boltzmannAverage, calibrate

# SMI, = glob_wildcards(abspath(join('input', '{id}.smi')))
# INCHI, = glob_wildcards(abspath(join('input', '{id}.inchi')))
# IDS = SMI + INCHI

# IDS.sort()

filepaths = glob.glob(abspath(join('input', '*.*')))
filepaths.sort()

d = {'mass': [], 'id': []}
for f in filepaths:
    ikey, ext = splitext(basename(f))
    d['mass'].append(next(pybel.readfile(ext[1:], f)).molwt)
    d['id'].append(ikey)

df = pd.DataFrame(d)
df.sort_values(by='mass', inplace=True)

IDS = df['id'].values

if 'stop' in config:
    IDS = IDS[config['start']:config['stop']]


rule all:
    input:
        expand(abspath(join('output', 'mobility', 'mobcal', 'calibrated_ccs', '{id}_{adduct}.tsv')),
               id=IDS, adduct=config['adducts'])
