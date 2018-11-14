from os.path import *
import numpy as np
import os
from isicle.utils import cycles, frames

# snakemake configuration
include: 'adducts.snakefile'
ruleorder: sander0 > sander


rule prepare:
    input:
        rules.generateAdduct.output.mol2
    output:
        mol2 = abspath(join('output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.input.mol2')),
        idx = abspath(join('output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.idx.npy')),
        content = abspath(join('output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.content.npy'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'md', 'antechamber', 'logs', '{id}_{adduct}.prepare.log'))
    benchmark:
        abspath(join('output', 'md', 'antechamber', 'benchmarks', '{id}_{adduct}.prepare.benchmark'))
    # group:
    #     'md'
    shell:
        # if charges come from DFT, use them (don't override with +1/-1/0)
        # also adjust antechamber flag if using DFT partial charges so it does not
        # assign
        'python -m isicle.scripts.md_helper {input} {output.mol2} [{wildcards.adduct}] --prepare &> {log}'


rule antechamber:
    input:
        mol2 = rules.prepare.output.mol2,
        idx = rules.prepare.output.idx
    output:
        mol2 = abspath(join('output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.output.mol2')),
        ac = abspath(join('output', 'md', 'antechamber', '{id}_{adduct}', 'ANTECHAMBER_AC.AC'))
    version:
        "antechamber | grep 'Welcome to antechamber' | awk '{print substr($4, 0, length($4) - 1)}'"
    log:
        abspath(join('output', 'md', 'antechamber', 'logs', '{id}_{adduct}.antechamber.log'))
    benchmark:
        abspath(join('output', 'md', 'antechamber', 'benchmarks', '{id}_{adduct}.antechamber.benchmark'))
    # group:
    #     'md'
    run:
        cwd = os.getcwd()
        os.chdir(abspath(join('output', 'md', 'antechamber', '%s_%s')) % (wildcards.id, wildcards.adduct))

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
        frcmod = abspath(join('output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.frcmod'))
    version:
        # using cpptraj as proxy for version
        "cpptraj --version | awk '{print substr($3, 2, length($3))}'"
    log:
        abspath(join('output', 'md', 'antechamber', 'logs', '{id}_{adduct}.parmchk2.log'))
    benchmark:
        abspath(join('output', 'md', 'antechamber', 'benchmarks', '{id}_{adduct}.parmchk2.benchmark'))
    # group:
    #     'md'
    run:
        cwd = os.getcwd()
        os.chdir(abspath(join('output', 'md', 'antechamber', '%s_%s')) % (wildcards.id, wildcards.adduct))

        shell('parmchk2 -i ANTECHAMBER_AC.AC -f ac -o {output.frcmod} &> {log}')

        os.chdir(cwd)


rule restore:
    input:
        mol2 = rules.antechamber.output.mol2,
        idx = rules.prepare.output.idx,
        content = rules.prepare.output.content
    output:
        mol2 = abspath(join('output', 'md', 'antechamber', '{id}_{adduct}', '{id}_{adduct}.mol2'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'md', 'antechamber', 'logs', '{id}_{adduct}.restore.log'))
    benchmark:
        abspath(join('output', 'md', 'antechamber', 'benchmarks', '{id}_{adduct}.restore.benchmark'))
    # group:
    #     'md'
    shell:
        'python -m isicle.scripts.md_helper {input.mol2} {output.mol2} [{wildcards.adduct}] --restore &> {log}'


rule tleapConfig:
    input:
        mol2 = rules.restore.output.mol2,
        frcmod = rules.parmchk2.output.frcmod
    output:
        abspath(join('output', 'md', 'tleap', '{id}_{adduct}.config'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'md', 'tleap', 'logs', '{id}_{adduct}.config.log'))
    benchmark:
        abspath(join('output', 'md', 'tleap', 'benchmarks', '{id}_{adduct}.config.benchmark'))
    # group:
    #     'md'
    shell:
        'python -m isicle.scripts.prepare_tleap {input.mol2} {input.frcmod} {output} &> {log}'


rule tleap:
    input:
        rules.tleapConfig.output
    output:
        prmtop = abspath(join('output', 'md', 'tleap', '{id}_{adduct}.top')),
        inpcrd = abspath(join('output', 'md', 'tleap', '{id}_{adduct}.crd'))
    version:
        # using cpptraj as proxy for version
        "cpptraj --version | awk '{print substr($3, 2, length($3))}'"
    log:
        abspath(join('output', 'md', 'tleap', 'logs', '{id}_{adduct}.log'))
    benchmark:
        abspath(join('output', 'md', 'tleap', 'benchmarks', '{id}_{adduct}.benchmark'))
    # group:
    #     'md'
    shell:
        'tleap -s -f {input} &> {log}'


rule sanderEMconfig:
    input:
        rules.restore.output.mol2
    output:
        abspath(join('output', 'md', 'em', '{id}_{adduct}.mdin'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'md', 'em', 'logs', '{id}_{adduct}.config.log'))
    benchmark:
        abspath(join('output', 'md', 'em', 'benchmarks', '{id}_{adduct}.config.benchmark'))
    # group:
    #     'md'
    shell:
        'python -m isicle.scripts.prepare_sander {input} {output} --em &> {log}'


rule sanderEM:
    input:
        prmtop = rules.tleap.output.prmtop,
        inpcrd = rules.tleap.output.inpcrd,
        config = rules.sanderEMconfig.output
    output:
        rst = abspath(join('output', 'md', 'em', '{id}_{adduct}.rst')),
        out = abspath(join('output', 'md', 'em', '{id}_{adduct}.out'))
    version:
        # using cpptraj as proxy for version
        "cpptraj --version | awk '{print substr($3, 2, length($3))}'"
    log:
        a = abspath(join('output', 'md', 'em', 'logs', '{id}_{adduct}.sander.log')),
        b = abspath(join('output', 'md', 'em', 'logs', '{id}_{adduct}.sander.log2'))
    benchmark:
        abspath(join('output', 'md', 'em', 'benchmarks', '{id}_{adduct}.sander.benchmark'))
    # group:
    #     'md'
    shell:
        'sander -O -i {input.config} -o {output.out} -c {input.inpcrd} -p {input.prmtop} -r {output.rst} \
         -inf {log.a} &> {log.b}'


rule sander0config:
    input:
        rules.restore.output.mol2
    output:
        abspath(join('output', 'md', 'anneal', 'cycle_000', '{id}_{adduct}.mdin'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'md', 'anneal', 'logs', '{id}_{adduct}_000.config.log'))
    benchmark:
        abspath(join('output', 'md', 'anneal', 'benchmarks', '{id}_{adduct}_000.config.benchmark'))
    # group:
    #     'md'
    shell:
        'python -m isicle.scripts.prepare_sander {input} {output} --iter0 &> {log}'


rule sander0:
    input:
        rst = rules.sanderEM.output.rst,
        prmtop = rules.tleap.output.prmtop,
        config = rules.sander0config.output
    output:
        rst = abspath(join('output', 'md', 'anneal', 'cycle_000', '{id}_{adduct}.rst')),
        crd = abspath(join('output', 'md', 'anneal', 'cycle_000', '{id}_{adduct}.crd')),
        out = abspath(join('output', 'md', 'anneal', 'cycle_000', '{id}_{adduct}.out'))
    version:
        # using cpptraj as proxy for version
        "cpptraj --version | awk '{print substr($3, 2, length($3))}'"
    log:
        a = abspath(join('output', 'md', 'anneal', 'logs', '{id}_{adduct}_000.sander.log')),
        b = abspath(join('output', 'md', 'anneal', 'logs', '{id}_{adduct}_000.sander.log2'))
    benchmark:
        abspath(join('output', 'md', 'anneal', 'benchmarks', '{id}_{adduct}_000.sander.benchmark'))
    # group:
    #     'md'
    shell:
        'sander -O -p {input.prmtop} -c {input.rst} -i {input.config} -o {output.out} -r {output.rst} \
         -x {output.crd} -inf {log.a} &> {log.b}'


rule sanderConfig:
    input:
        rules.restore.output.mol2
    output:
        abspath(join('output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.mdin'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'md', 'anneal', 'logs', '{id}_{adduct}_{cycle}.config.log'))
    benchmark:
        abspath(join('output', 'md', 'anneal', 'benchmarks', '{id}_{adduct}_{cycle}.config.benchmark'))
    # group:
    #     'md'
    shell:
        'python -m isicle.scripts.prepare_sander {input} {output} --sa &> {log}'


rule sander:
    input:
        # s0 required to disambiguate, but not used
        rst0 = rules.sander0.output.rst,
        rst = lambda wildcards: abspath(join('output', 'md', 'anneal', 'cycle_%03d', '%s_%s.rst')) % \
                                (int(wildcards.cycle) - 1, wildcards.id, wildcards.adduct),
        prmtop = rules.tleap.output.prmtop,
        config = rules.sanderConfig.output
    output:
        rst = abspath(join('output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.rst')),
        crd = abspath(join('output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.crd')),
        out = abspath(join('output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.out'))
    version:
        # using cpptraj as proxy for version
        "cpptraj --version | awk '{print substr($3, 2, length($3))}'"
    log:
        a = abspath(join('output', 'md', 'anneal', 'logs', '{id}_{adduct}_{cycle}.sander.log')),
        b = abspath(join('output', 'md', 'anneal', 'logs', '{id}_{adduct}_{cycle}.sander.log2'))
    benchmark:
        abspath(join('output', 'md', 'anneal', 'benchmarks', '{id}_{adduct}_{cycle}.sander.benchmark'))
    # group:
    #     'md'
    shell:
        'sander -O -p {input.prmtop} -c {input.rst} -i {input.config} -o {output.out} -r {output.rst} \
         -x {output.crd} -inf {log.a} &> {log.b}'


rule selectFrames:
    input:
        out = abspath(join('output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.out')),
        crd = abspath(join('output', 'md', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.crd'))
    output:
        expand(abspath(join('output', 'md', 'extracted', '{{id}}_{{adduct}}_{{cycle}}_{frame}.trajin')),
               frame=frames(config['amber']['nframes']))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'md', 'extracted', 'logs', '{id}_{adduct}_{cycle}.select.log'))
    benchmark:
        abspath(join('output', 'md', 'extracted', 'benchmarks', '{id}_{adduct}_{cycle}.select.benchmark'))
    # group:
    #     'md'
    shell:
        'python -m isicle.scripts.select_frames {input.out} {input.crd} {output} \
         --nframes {config[amber][nframes]} --low {config[amber][low]} --high {config[amber][high]} &> {log}'


rule extractFrames:
    input:
        prmtop = rules.tleap.output.prmtop,
        trajin = abspath(join('output', 'md', 'extracted', '{id}_{adduct}_{cycle}_{frame}.trajin'))
    output:
        abspath(join('output', 'md', 'extracted', '{id}_{adduct}_{cycle}_{frame}.mol2'))
    version:
        "cpptraj --version | awk '{print substr($3, 2, length($3))}'"
    log:
        abspath(join('output', 'md', 'extracted', 'logs', '{id}_{adduct}_{cycle}_{frame}.extract.log'))
    benchmark:
        abspath(join('output', 'md', 'extracted', 'benchmarks', '{id}_{adduct}_{cycle}_{frame}.extract.benchmark'))
    # group:
    #     'md'
    shell:
        'cpptraj {input.prmtop} {input.trajin} &> {log}'


rule convert:
    input:
        mol2a = rules.extractFrames.output,
        mol2b = rules.prepare.input
    output:
        abspath(join('output', 'md', 'converted', '{id}_{adduct}_{cycle}_{frame}.xyz'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'md', 'converted', 'logs', '{id}_{adduct}_{cycle}_{frame}.log'))
    benchmark:
        abspath(join('output', 'md', 'converted', 'benchmarks', '{id}_{adduct}_{cycle}_{frame}.benchmark'))
    # group:
    #     'md'
    shell:
        'python -m isicle.scripts.standardize_mol2 {input.mol2a} {input.mol2b} {output} &> {log}'


rule calculate_rmsd:
    input:
        ref = abspath(join('output', 'md', 'converted', '{id}_{adduct}_{cycle}_{frame}.xyz')),
        xyzs = expand(abspath(join('output', 'md', 'converted', '{{id}}_{{adduct}}_{{cycle}}_{frame}.xyz')),
                      frame=frames(config['amber']['nframes']))
    output:
        rmsd = abspath(join('output', 'md', 'rmsd', '{id}_{adduct}_{cycle}_{frame}.rmsd'))
    version:
        'isicle --version'
    log:
        abspath(join('output', 'md', 'rmsd', 'logs', '{id}_{adduct}_{cycle}_{frame}.log'))
    benchmark:
        abspath(join('output', 'md', 'rmsd', 'benchmarks', '{id}_{adduct}_{cycle}_{frame}.benchmark'))
    # group:
    #     'md'
    shell:
        'python -m isicle.scripts.rmsd {input.ref} {output.rmsd} {input.xyzs} &> {log}'


rule downselect:
    input:
        xyz = expand(abspath(join('output', 'md', 'converted', '{{id}}_{{adduct}}_{{cycle}}_{frame}.xyz')),
                     frame=frames(config['amber']['nframes'])),
        rmsd = expand(abspath(join('output', 'md', 'rmsd', '{{id}}_{{adduct}}_{{cycle}}_{frame}.rmsd')),
                      frame=frames(config['amber']['nframes']))
    output:
        expand(abspath(join('output', 'md', 'downselected', '{{id}}_{{adduct}}_{{cycle}}_{selected}.xyz')),
               selected=['s', 'd1', 'd2'])
    version:
        'isicle --version'
    log:
        abspath(join('output', 'md', 'downselected', 'logs', '{id}_{adduct}_{cycle}.log'))
    benchmark:
        abspath(join('output', 'md', 'downselected', 'benchmarks', '{id}_{adduct}_{cycle}.benchmark'))
    # group:
    #     'md'
    run:
        outdir = abspath(join('output', 'md', 'downselected'))
        shell('python -m isicle.scripts.downselect {outdir} --infiles {input.xyz} --rfiles {input.rmsd} &> {log}')
