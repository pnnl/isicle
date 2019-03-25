from os.path import *
from isicle.utils import cycles
import shutil
from pkg_resources import resource_filename
import os

# snakemake configuration
include: 'shielding.snakefile'
# localrules: all, inchi2smiles, canonicalize, desalt, neutralize, tautomerize,
#             calculateFormula, calculateMass, generateGeometry, generateAdduct,
#             prepare, restore, tleapConfig, sanderEMconfig, sander0config, sanderConfig,
#             selectFrames, extractFrames, convert, calculate_rmsd, downselect, copyOver,
#             createShieldingConfig, parseShielding, combine, boltzmannAverage, shifts

SMI, = glob_wildcards(abspath(join('input', '{id}.smi')))
INCHI, = glob_wildcards(abspath(join('input', '{id}.inchi')))
IDS = SMI + INCHI

IDS.sort()
if 'stop' in config:
    IDS = IDS[config['start']:config['stop']]

rule all:
    input:
        expand(abspath(join('output', 'shifts', '{id}.tsv')), id=IDS)
