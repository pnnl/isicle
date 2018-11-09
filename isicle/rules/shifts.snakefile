from os.path import *
from isicle.utils import cycles
import shutil
from pkg_resources import resource_filename
import os

# snakemake configuration
include: 'shielding.snakefile'

SMI, = glob_wildcards(abspath(join('input', '{id}.smi')))
INCHI, = glob_wildcards(abspath(join('input', '{id}.inchi')))
IDS = SMI + INCHI

# copy reference molecule
if config['nwchem']['reference'] in ['TMS', 'DSS']:
    if not exists(abspath('input')):
        os.mkdir(abspath('input'))
    if not exists(abspath(join('input', config['nwchem']['reference'] + '.smi'))):
        shutil.copy2(resource_filename('isicle', join('resources', 'nwchem', config['nwchem']['reference'] + '.smi')),
                     abspath(join('input', config['nwchem']['reference'] + '.smi')))
        IDS.append(config['nwchem']['reference'])

# check if supplied reference is in input
elif config['nwchem']['reference'] not in IDS:
    raise Exception('Select TMS or DSS as reference, or ensure %s is in the input folder.' % config['nwchem']['reference'])

IDS.sort()
if 'stop' in config:
    IDS = IDS[config['start']:config['stop']]

# add reference back in
if config['nwchem']['reference'] not in IDS:
    IDS.append(config['nwchem']['reference'])


rule all:
    input:
        expand(abspath(join('output', 'shifts', '{id}.tsv')), id=IDS)
