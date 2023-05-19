from os.path import *

# snakemake configuration
include: 'mobility_alt.snakefile'

SMI, = glob_wildcards(abspath(join('input', '{id}.smi')))
INCHI, = glob_wildcards(abspath(join('input', '{id}.inchi')))
IDS = SMI + INCHI

IDS.sort()

if 'stop' in config:
    IDS = IDS[config['start']:config['stop']]


rule all:
    input:
        expand(abspath(join('combined_ccs','{id}_{adduct}.He.txt')),
               id=IDS, adduct=config['adducts']),
        expand(abspath(join('combined_ccs','{id}_{adduct}.N2.txt')),
               id=IDS, adduct=config['adducts'])

def aggregate_adducts_He(wildcards):
    checkpoint_output = checkpoints.generateAdducts.get(**wildcards).output[0]
    return expand(abspath(join('output', 'mobility', 'impact', 'ccs', '{id}_{adduct}_{addID}.He.ccs')),
        id=wildcards.id,
        adduct=wildcards.adduct,
        addID=glob_wildcards(join(checkpoint_output,'{addID}.charge')).addID)

def aggregate_adducts_N2(wildcards):
    checkpoint_output = checkpoints.generateAdducts.get(**wildcards).output[0]
    return expand(abspath(join('output', 'mobility', 'impact', 'ccs', '{id}_{adduct}_{addID}.N2.ccs')),
        id=wildcards.id,
        adduct=wildcards.adduct,
        addID=glob_wildcards(join(checkpoint_output,'{addID}.charge')).addID)

rule collectOutput:
    input:
        heccs= aggregate_adducts_He,
        n2ccs= aggregate_adducts_N2
    output:
        heccs= join('combined_ccs','{id}_{adduct}.He.txt')
        n2ccs= join('combined_ccs','{id}_{adduct}.N2.txt')
    shell:
        '''
        cat {input.heccs} > {output.heccs}
        cat {input.n2ccs} > {output.n2ccs}
        '''