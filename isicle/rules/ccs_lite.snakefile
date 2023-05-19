from os.path import *

# snakemake configuration
include: 'adducts.snakefile'
include: 'mobility_alt.snakefile'
localrules: all, inchi2smiles, canonicalize, desalt, neutralize, calculateMass, generateGeometry,
            touchAdducts, collectOutput, impact, postprocess

SMI, = glob_wildcards(abspath(join('input', '{id}.smi')))
INCHI, = glob_wildcards(abspath(join('input', '{id}.inchi')))
IDS = SMI + INCHI

IDS.sort()

if 'stop' in config:
    IDS = IDS[config['start']:config['stop']]

rule all:
    input:
        expand(join('output', 'combined_ccs','{id}_{adduct}.txt'),
               id=IDS, adduct=config['adducts'])

def aggregate_adducts(wildcards):
    checkpoint_output = checkpoints.generateAdducts.get(**wildcards).output[0]
    heccs= expand(abspath(join('output', 'mobility', 'impact', 'ccs', '{id}_{adduct}_{addID}.He.ccs')),
        id=wildcards.id,
        adduct=wildcards.adduct,
        addID=glob_wildcards(join(checkpoint_output,'{addID}.charge')).addID)
    n2ccs= expand(abspath(join('output', 'mobility', 'impact', 'ccs', '{id}_{adduct}_{addID}.N2.ccs')),
        id=wildcards.id,
        adduct=wildcards.adduct,
        addID=glob_wildcards(join(checkpoint_output,'{addID}.charge')).addID)
    return heccs+n2ccs

rule collectOutput:
    input:
        aggregate_adducts
    output:
        join('output', 'combined_ccs','{id}_{adduct}.txt')
    shell:
        '''
        echo {input} >> {output}
        cat {input} >> {output}
        '''