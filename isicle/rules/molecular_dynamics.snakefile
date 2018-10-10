from os.path import *
import numpy as np
import os
from generators import *

# snakemake configuration
include: 'adducts.snakefile'
ruleorder: sander0 > sander


rule prepare:
    input:
        rules.generateAdduct.output.mol2
    output:
        mol2 = join(config['path'], 'output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.input.mol2'),
        idx = join(config['path'], 'output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.idx.npy'),
        content = join(config['path'], 'output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.content.npy')
    log:
        join(config['path'], 'output', 'md', 'antechamber', 'logs', '{id}_{adduct}.prepare.log')
    benchmark:
        join(config['path'], 'output', 'md', 'antechamber', 'benchmarks', '{id}_{adduct}.prepare.benchmark')
    # group:
    #     'md'
    shell:
        # if charges come from DFT, use them (don't override with +1/-1/0)
        # also adjust antechamber flag if using DFT partial charges so it does not
        # assign
        'python isicle/md_helper.py {input} {output.mol2} [{wildcards.adduct}] --prepare &> {log}'


rule antechamber:
    input:
        mol2 = rules.prepare.output.mol2,
        idx = rules.prepare.output.idx
    output:
        mol2 = join(config['path'], 'output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.output.mol2'),
        ac = join(config['path'], 'output', 'md', 'antechamber', '{id}_{adduct}', 'ANTECHAMBER_AC.AC')
    log:
        join(config['path'], 'output', 'md', 'antechamber', 'logs', '{id}_{adduct}.antechamber.log')
    benchmark:
        join(config['path'], 'output', 'md', 'antechamber', 'benchmarks', '{id}_{adduct}.antechamber.benchmark')
    # group:
    #     'md'
    run:
        cwd = os.getcwd()
        os.chdir(join(config['path'], 'output', 'md', 'antechamber', '%s_%s') % (wildcards.id, wildcards.adduct))

        charge = config['charges'][wildcards.adduct]
        if wildcards.adduct == '+Na':
            idx = np.load(input.idx)
            charge -= len(idx)

        shell('antechamber -i {input.mol2} -fi mol2 -o {output.mol2} -fo mol2 -c bcc -s -du -nc %.4f &> {log}' % charge)

        os.chdir(cwd)


rule parmchk2:
    input:
        rules.antechamber.output.ac
    output:
        frcmod = join(config['path'], 'output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.frcmod')
    log:
        join(config['path'], 'output', 'md', 'antechamber', 'logs', '{id}_{adduct}.parmchk2.log')
    benchmark:
        join(config['path'], 'output', 'md', 'antechamber', 'benchmarks', '{id}_{adduct}.parmchk2.benchmark')
    # group:
    #     'md'
    run:
        cwd = os.getcwd()
        os.chdir(join(config['path'], 'output', 'md', 'antechamber', '%s_%s') % (wildcards.id, wildcards.adduct))

        shell('parmchk2 -i ANTECHAMBER_AC.AC -f ac -o {output.frcmod} &> {log}')

        os.chdir(cwd)


rule restore:
    input:
        mol2 = rules.antechamber.output.mol2,
        idx = rules.prepare.output.idx,
        content = rules.prepare.output.content
    output:
        mol2 = join(config['path'], 'output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.mol2')
    log:
        join(config['path'], 'output', 'md', 'antechamber', 'logs', '{id}_{adduct}.restore.log')
    benchmark:
        join(config['path'], 'output', 'md', 'antechamber', 'benchmarks', '{id}_{adduct}.restore.benchmark')
    # group:
    #     'md'
    shell:
        'python isicle/md_helper.py {input.mol2} {output.mol2} [{wildcards.adduct}] --restore &> {log}'


rule tleap:
    input:
        mol2 = rules.restore.output.mol2,
        frcmod = rules.parmchk2.output.frcmod
    output:
        config = join(config['path'], 'output', 'md', 'tleap', '{id}_{adduct}.config'),
        prmtop = join(config['path'], 'output', 'md', 'tleap', '{id}_{adduct}.top'),
        inpcrd = join(config['path'], 'output', 'md', 'tleap', '{id}_{adduct}.crd')
    log:
        join(config['path'], 'output', 'md', 'tleap', 'logs', '{id}_{adduct}.meta.log')
    benchmark:
        join(config['path'], 'output', 'md', 'tleap', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'md'
    run:
        shell('python isicle/prepare_tleap.py {input.mol2} {input.frcmod} {output.config} &> {log}')
        shell('tleap -s -f {output.config} &> {log}')


rule sanderEM:
    input:
        mol2 = rules.restore.output.mol2,
        prmtop = rules.tleap.output.prmtop,
        inpcrd = rules.tleap.output.inpcrd
    output:
        config = join(config['path'], 'output', 'md', 'em', '{id}_{adduct}.mdin'),
        rst = join(config['path'], 'output', 'md', 'em', '{id}_{adduct}.rst'),
        out = join(config['path'], 'output', 'md', 'em', '{id}_{adduct}.out')
    log:
        join(config['path'], 'output', 'md', 'em', 'logs', '{id}_{adduct}.log')
    benchmark:
        join(config['path'], 'output', 'md', 'em', 'benchmarks', '{id}_{adduct}.benchmark')
    # group:
    #     'md'
    run:
        shell('python isicle/prepare_sander.py {input.mol2} {output.config} --em &> {log}')
        shell('sander -O -i {output.config} -o {output.out} -c {input.inpcrd} -p {input.prmtop} -r {output.rst} -inf {log}')


rule sander0:
    input:
        mol2 = rules.restore.output.mol2,
        rst = rules.sanderEM.output.rst,
        prmtop = rules.tleap.output.prmtop
    output:
        config = join(config['path'], 'output', 'md', 'anneal', 'cycle_000', '{id}_{adduct}.mdin'),
        rst = join(config['path'], 'output', 'md', 'anneal', 'cycle_000', '{id}_{adduct}.rst'),
        crd = join(config['path'], 'output', 'md', 'anneal', 'cycle_000', '{id}_{adduct}.crd'),
        out = join(config['path'], 'output', 'md', 'anneal', 'cycle_000', '{id}_{adduct}.out')
    log:
        join(config['path'], 'output', 'md', 'anneal', 'logs', '{id}_{adduct}_000.log')
    benchmark:
        join(config['path'], 'output', 'md', 'anneal', 'benchmarks', '{id}_{adduct}_000.benchmark')
    # group:
    #     'md'
    run:
        shell('python isicle/prepare_sander.py {input.mol2} {output.config} --iter0 &> {log}')
        shell('sander -O -p {input.prmtop} -c {input.rst} -i {output.config} -o {output.out} -r {output.rst} -x {output.crd} -inf {log}')


rule sander:
    input:
        mol2 = rules.restore.output.mol2,
        # s0 required to disambiguate, but not used
        rst0 = rules.sander0.output.rst,
        rst = lambda wildcards: join(config['path'], 'output', 'md', 'anneal', 'cycle_%03d', '%s_%s.rst') % \
                                    (int(wildcards.cycle) - 1, wildcards.id, wildcards.adduct),
        prmtop = rules.tleap.output.prmtop
    output:
        config = join(config['path'], 'output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.mdin'),
        rst = join(config['path'], 'output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.rst'),
        crd = join(config['path'], 'output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.crd'),
        out = join(config['path'], 'output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.out')
    log:
        join(config['path'], 'output', 'md', 'anneal', 'logs', '{id}_{adduct}_{cycle}.log')
    benchmark:
        join(config['path'], 'output', 'md', 'anneal', 'benchmarks', '{id}_{adduct}_{cycle}.benchmark')
    # group:
    #     'md'
    run:
        shell('python isicle/prepare_sander.py {input.mol2} {output.config} --sa &> {log}')
        shell('sander -O -p {input.prmtop} -c {input.rst} -i {output.config} -o {output.out} -r {output.rst} -x {output.crd} -inf {log}')


rule selectFrames:
    input:
        out = join(config['path'], 'output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.out'),
        crd = join(config['path'], 'output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.crd')
    output:
        expand(join(config['path'], 'output', 'md', 'extracted', '{{id}}_{{adduct}}_{{cycle}}_{frame}.trajin'),
               frame=frames(config['amber']['nframes']))
    log:
        join(config['path'], 'output', 'md', 'extracted', 'logs', '{id}_{adduct}_{cycle}.select.log')
    benchmark:
        join(config['path'], 'output', 'md', 'extracted', 'benchmarks', '{id}_{adduct}_{cycle}.select.benchmark')
    # group:
    #     'md'
    shell:
        'python isicle/select_frames.py {input.out} {input.crd} {output} \
         --nframes {config[amber][nframes]} --low {config[amber][low]} --high {config[amber][high]} &> {log}'


rule extractFrames:
    input:
        prmtop = rules.tleap.output.prmtop,
        trajin = join(config['path'], 'output', 'md', 'extracted', '{id}_{adduct}_{cycle}_{frame}.trajin')
    output:
        join(config['path'], 'output', 'md', 'extracted', '{id}_{adduct}_{cycle}_{frame}.mol2')
    log:
        join(config['path'], 'output', 'md', 'extracted', 'logs', '{id}_{adduct}_{cycle}_{frame}.extract.log')
    benchmark:
        join(config['path'], 'output', 'md', 'extracted', 'benchmarks', '{id}_{adduct}_{cycle}_{frame}.extract.benchmark')
    # group:
    #     'md'
    shell:
        'cpptraj {input.prmtop} {input.trajin} &> {log}'


rule convert:
    input:
        mol2a = rules.extractFrames.output,
        mol2b = rules.prepare.input
    output:
        join(config['path'], 'output', 'md', 'converted', '{id}_{adduct}_{cycle}_{frame}.xyz')
    log:
        join(config['path'], 'output', 'md', 'converted', 'logs', '{id}_{adduct}_{cycle}_{frame}.log')
    benchmark:
        join(config['path'], 'output', 'md', 'converted', 'benchmarks', '{id}_{adduct}_{cycle}_{frame}.benchmark')
    # group:
    #     'md'
    shell:
        'python isicle/standardize_mol2.py {input.mol2a} {input.mol2b} {output} &> {log}'


rule calculate_rmsd:
    input:
        ref = join(config['path'], 'output', 'md', 'converted', '{id}_{adduct}_{cycle}_{frame}.xyz'),
        xyzs = expand(join(config['path'], 'output', 'md', 'converted', '{{id}}_{{adduct}}_{{cycle}}_{frame}.xyz'),
                      frame=frames(config['amber']['nframes']))
    output:
        rmsd = join(config['path'], 'output', 'md', 'rmsd', '{id}_{adduct}_{cycle}_{frame}.rmsd')
    log:
        join(config['path'], 'output', 'md', 'rmsd', 'logs', '{id}_{adduct}_{cycle}_{frame}.log')
    benchmark:
        join(config['path'], 'output', 'md', 'rmsd', 'benchmarks', '{id}_{adduct}_{cycle}_{frame}.benchmark')
    # group:
    #     'md'
    shell:
        'python isicle/rmsd.py {input.ref} {output.rmsd} {input.xyzs} &> {log}'


rule downselect:
    input:
        xyz = expand(join(config['path'], 'output', 'md', 'converted', '{{id}}_{{adduct}}_{{cycle}}_{frame}.xyz'),
                     frame=frames(config['amber']['nframes'])),
        rmsd = expand(join(config['path'], 'output', 'md', 'rmsd', '{{id}}_{{adduct}}_{{cycle}}_{frame}.rmsd'),
                      frame=frames(config['amber']['nframes']))
    output:
        expand(join(config['path'], 'output', 'md', 'downselected', '{{id}}_{{adduct}}_{{cycle}}_{selected}.xyz'),
               selected=['s', 'd1', 'd2'])
    log:
        join(config['path'], 'output', 'md', 'downselected', 'logs', '{id}_{adduct}_{cycle}.log')
    benchmark:
        join(config['path'], 'output', 'md', 'downselected', 'benchmarks', '{id}_{adduct}_{cycle}.benchmark')
    # group:
    #     'md'
    run:
        shell('python isicle/downselect.py %s --infiles {input.xyz} --rfiles {input.rmsd} &> {log}' %
              join(config['path'], 'output', 'md', 'downselected'))
