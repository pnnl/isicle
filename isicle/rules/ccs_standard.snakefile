from os.path import *
from openbabel import pybel
import pandas as pd
import glob

# snakemake configuration
include: 'adducts.snakefile'
include: 'mobility.snakefile'
localrules: all, inchi2smiles, canonicalize, desalt, neutralize, calculateMass, generateGeometry,
            touchAdducts, prepare, restore, tleapConfig, sanderEMconfig, sander0config, sanderConfig,
            selectFrames, extractFrames, convert, calculate_rmsd, downselect, copyOver,
            createDFTConfig, parseDFT, parseMobcal, boltzmannAverage, calibrate, collectOutput

SMI, = glob_wildcards(abspath(join('input', '{id}.smi')))
INCHI, = glob_wildcards(abspath(join('input', '{id}.inchi')))
IDS = SMI + INCHI

IDS.sort()

print(IDS)


rule all:
    input:
        expand(join('output', 'combined_ccs','{id}_{adduct}.txt'),
               id=IDS, adduct=config['adducts'])

def aggregate_adducts(wildcards):
    checkpoint_output = checkpoints.generateAdducts.get(**wildcards).output[0]
    return expand(abspath(join('output', 'mobility', 'mobcal', 'calibrated_ccs', '{id}_{adduct}_{addID}.tsv')),
        id=wildcards.id,
        adduct=wildcards.adduct,
        addID=glob_wildcards(join(checkpoint_output,'{addID}.charge')).addID)

rule collectOutput:
    input:
        aggregate_adducts
    output:
        join('output', 'combined_ccs','{id}_{adduct}.txt')
    shell:
        '''
        cat {input} > {output}
        '''