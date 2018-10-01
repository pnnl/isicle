from os.path import *
from resources.utils import *
import shutil
import subprocess
import os
from string import Template
import glob
import numpy as np

# snakemake configuration
configfile: 'config.yaml'
include: 'adducts.snakefile'
ruleorder: sander0 > sander


rule prepare:
    input:
        mol2 = rules.generateAdducts.output.mol2
    output:
        mol2 = join(config['path'], 'output', 'antechamber', 'tmp', '{id}_{adduct}', '{id}_{adduct}.input.mol2'),
        idx = join(config['path'], 'output', 'antechamber', 'tmp', '{id}_{adduct}', '{id}_{adduct}.idx.npy'),
        content = join(config['path'], 'output', 'antechamber', 'tmp', '{id}_{adduct}', '{id}_{adduct}.content.npy')
    group:
        'md'
    run:
        # if charges come from DFT, use them (don't override with +1/-1/0)
        # also adjust antechamber flag if using DFT partial charges so it does not
        # assign
        mol = read_mol(input.mol2)

        if wildcards.adduct == '+Na':
            idx, content = pop_atom(input.mol2, output.mol2, atom='Na')
        else:
            shutil.copy2(input.mol2, output.mol2)
            idx = None
            content = None

        np.save(output.idx, idx)
        np.save(output.content, content)

rule antechamber:
    input:
        mol2 = rules.prepare.output.mol2,
        idx = rules.prepare.output.idx
    output:
        mol2 = join(config['path'], 'output', 'antechamber', 'tmp', '{id}_{adduct}', '{id}_{adduct}.output.mol2'),
        ac = join(config['path'], 'output', 'antechamber', 'tmp', '{id}_{adduct}', 'ANTECHAMBER_AC.AC')
    log:
        join(config['path'], 'output', 'antechamber', 'logs', '{id}_{adduct}.antechamber.log'),
    group:
        'md'
    run:
        cwd = os.getcwd()
        os.chdir(join(config['path'], 'output', 'antechamber', 'tmp', '%s_%s' % (wildcards.id, wildcards.adduct)))

        charge = config['charges'][wildcards.adduct]
        if wildcards.adduct == '+Na':
            idx = np.load(input.idx)
            charge -= len(idx)

        shell('antechamber -i {input.mol2} -fi mol2 -o {output.mol2} -fo mol2 -c bcc -s -du -nc %.4f &> {log}' % charge)
        os.chdir(cwd)

rule parmchk2:
    input:
        ac = rules.antechamber.output.mol2
    output:
        frcmod = join(config['path'], 'output', 'antechamber', 'frcmod', '{id}_{adduct}.frcmod')
    log:
        join(config['path'], 'output', 'antechamber', 'logs', '{id}_{adduct}.parmchk2.log')
    group:
        'md'
    run:
        cwd = os.getcwd()
        os.chdir(join(config['path'], 'output', 'antechamber', 'tmp', '%s_%s' % (wildcards.id, wildcards.adduct)))

        shell('parmchk2 -i ANTECHAMBER_AC.AC -f ac -o {output.frcmod} &> {log}')

        os.chdir(cwd)

rule restore:
    input:
        mol2 = rules.antechamber.output.mol2,
        idx = rules.prepare.output.idx,
        content = rules.prepare.output.content
    output:
        mol2 = join(config['path'], 'output', 'antechamber', 'mol2', '{id}_{adduct}.mol2')
    group:
        'md'
    run:
        if wildcards.adduct == '+Na':
            idx = np.load(input.idx)
            content = np.load(input.content)
            push_atom(input.mol2, output.mol2, idx, content)
        else:
            shutil.copy2(input.mol2, output.mol2)

rule tleapConfig:
    input:
        mol2 = rules.antechamber.output.mol2,
        frcmod = rules.parmchk2.output.frcmod
    output:
        config = join(config['path'], 'output', 'tleap', 'config', '{id}_{adduct}.config')
    group:
        'md'
    run:
        with open('resources/amber/tleap.template', 'r') as f:
            t = Template(f.read())

        d = {'frcmod': input.frcmod,
             'mol2': input.mol2,
             'prmtop': join(config['path'], 'output', 'tleap', 'prmtop', '%s_%s.top' % (wildcards.id, wildcards.adduct)),
             'inpcrd': join(config['path'], 'output', 'tleap', 'inpcrd', '%s_%s.crd' % (wildcards.id, wildcards.adduct))}

        with open(output.config, 'w') as f:
            f.write(t.substitute(d))

rule tleap:
    input:
        config = rules.tleapConfig.output.config
    output:
        prmtop = join(config['path'], 'output', 'tleap', 'prmtop', '{id}_{adduct}.top'),
        inpcrd = join(config['path'], 'output', 'tleap', 'inpcrd', '{id}_{adduct}.crd')
    log:
        join(config['path'], 'output', 'tleap', 'logs', '{id}_{adduct}.log')
    group:
        'md'
    shell:
        'tleap -s -f {input.config} logFile {log}'

rule sanderEMConfig:
    input:
        mol2 = rules.antechamber.output.mol2
    output:
        config = join(config['path'], 'output', 'sander', 'em', '{id}_{adduct}.mdin')
    group:
        'md'
    run:
        with open('resources/amber/sander_em.template', 'r') as f:
                t = Template(f.read())

        d = {'mol2': input.mol2,
             'imin': 1,
             'maxcyc': 500,
             'ncyc': 250,
             'nscm': 1,
             'ntb': 0,
             'igb': 0,
             'cut': 999,
             'rgbmax': 999}

        with open(output.config, 'w') as f:
            f.write(t.substitute(d))

rule sanderEM:
    input:
        config = rules.sanderEMConfig.output.config,
        mol2 = rules.restore.output.mol2,
        prmtop = rules.tleap.output.prmtop,
        inpcrd = rules.tleap.output.inpcrd
    output:
        rst = join(config['path'], 'output', 'sander', 'em', '{id}_{adduct}.rst'),
        out = join(config['path'], 'output', 'sander', 'em', '{id}_{adduct}.out')
    group:
        'md'
    shell:
        'sander -O -i {input.config} -o {output.out} -c {input.inpcrd} -p {input.prmtop} -r {output.rst}'

rule sander0:
    input:
        mol2 = rules.restore.output.mol2,
        rst = rules.sanderEM.output.rst,
        prmtop = rules.tleap.output.prmtop
    output:
        config = join(config['path'], 'output', 'sander', 'anneal', 'cycle_000', '{id}_{adduct}.mdin'),
        rst = join(config['path'], 'output', 'sander', 'anneal', 'cycle_000', '{id}_{adduct}.rst'),
        crd = join(config['path'], 'output', 'sander', 'anneal', 'cycle_000', '{id}_{adduct}.crd'),
        out = join(config['path'], 'output', 'sander', 'anneal', 'cycle_000', '{id}_{adduct}.out')
    group:
        'md'
    run:
        with open('resources/amber/sander_md0.template', 'r') as f:
            t = Template(f.read())

        d = {'mol2': input.mol2,
             'imin': 0,
             'ntb': 0,
             'ntf': 2,
             'ntc': 2,
             'ntx': 1,
             'igb': 0,
             'ntpr': 100,
             'ntwx': 100,
             'ntt': 3,
             'nscm': 1,
             'gamma_ln': 1.0,
             'tempi': 0.0,
             'temp0': 300.0,
             'nstlim': 100000,
             'dt': 0.0005,
             'cut': 999,
             'ig': -1}

        with open(output.config, 'w') as f:
            f.write(t.substitute(d))

        shell('sander -O -p {input.prmtop} -c {input.rst} -i {output.config} -o {output.out} -r {output.rst} -x {output.crd}')

rule sander:
    input:
        mol2 = rules.restore.output.mol2,
        # s0 required to disambiguate, but not used
        rst0 = rules.sander0.output.rst,
        rst = lambda wildcards: join(config['path'], 'output', 'sander', 'anneal', 'cycle_%03d', '%s_%s.rst') %
                                    (int(wildcards.cycle) - 1, wildcards.id, wildcards.adduct),
        prmtop = rules.tleap.output.prmtop
    output:
        config = join(config['path'], 'output', 'sander', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.mdin'),
        rst = join(config['path'], 'output', 'sander', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.rst'),
        crd = join(config['path'], 'output', 'sander', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.crd'),
        out = join(config['path'], 'output', 'sander', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.out')
    group:
        'md'
    run:
        with open('resources/amber/sander_anneal.template', 'r') as f:
            t = Template(f.read())

        d = {'mol2': input.mol2}

        with open(output.config, 'w') as f:
            f.write(t.substitute(d))

        shell('sander -O -p {input.prmtop} -c {input.rst} -i {output.config} -o {output.out} -r {output.rst} -x {output.crd}')

rule extractFrames:
    input:
        prmtop = rules.tleap.output.prmtop,
        out = join(config['path'], 'output', 'sander', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.out'),
        crd = join(config['path'], 'output', 'sander', 'anneal', 'cycle_{cycle}', '{id}_{adduct}.crd')
    output:
        trajin = join(config['path'], 'output', 'sander', 'extracted', 'trajin', '{id}_{adduct}_{cycle}_{frame}.trajin'),
        mol2 = join(config['path'], 'output', 'sander', 'extracted', 'mol2', '{id}_{adduct}_{cycle}_{frame}.mol2')
    log:
        join(config['path'], 'output', 'sander', 'extracted', 'logs', '{id}_{adduct}_{cycle}_{frame}.log')
    group:
        'md'
    run:
        frame = select_frames(input.out,
                              frames=config['amber']['nframes'],
                              low=config['amber']['low'],
                              high=config['amber']['high'])[int(wildcards.frame)]

        shell('echo "trajin {input.crd} %s %s" > {output.trajin}' % (frame, frame))
        shell('echo "trajout {output.mol2} mol2" >> {output.trajin}')
        shell('cpptraj {input.prmtop} {output.trajin} > {log}')

rule convert:
    input:
        mol2a = join(config['path'], 'output', 'sander', 'extracted', 'mol2', '{id}_{adduct}_{cycle}_{frame}.mol2'),
        mol2b = rules.antechamber.input.mol2
    output:
        xyz = join(config['path'], 'output', 'sander', 'extracted', 'xyz', '{id}_{adduct}_{cycle}_{frame}.xyz')
    group:
        'md'
    run:
        standardizeMol2(input.mol2a, input.mol2b, output.xyz)

rule calculate_rmsd:
    input:
        xyz = rules.convert.output.xyz
    output:
        rmsd = join(config['path'], 'output', 'sander', 'extracted', 'rmsd', '{id}_{adduct}_{cycle}_{frame}.rmsd')
    group:
        'md'
    run:
        mols = glob.glob(join(config['path'], 'output', 'sander', 'extracted', 'xyz', '%s_%s_%s_*.xyz' %
                              (wildcards.id, wildcards.adduct, wildcards.cycle)))
        total = 0
        for mol in mols:
            total += rmsd(input.xyz, mol)

        with open(output.rmsd, 'w') as f:
            f.write(str(total) + '\n')

rule downselect:
    input:
        xyz = expand(join(config['path'], 'output', 'sander', 'extracted', 'xyz', '{{id}}_{{adduct}}_{{cycle}}_{frame}.xyz'),
                     frame=frames(config['amber']['nframes'])),
        rmsd = expand(join(config['path'], 'output', 'sander', 'extracted', 'rmsd', '{{id}}_{{adduct}}_{{cycle}}_{frame}.rmsd'),
                      frame=frames(config['amber']['nframes']))
    output:
        selected = expand(join(config['path'], 'output', 'nwchem', '{{id}}_{{adduct}}_{{cycle}}_{selected}.xyz'),
                          selected=['s', 'd1', 'd2'])
    group:
        'md'
    run:
        vals = []
        for f in input.rmsd:
            vals.append(float(read_string(f)))

        vals = np.array(vals)
        idx = np.argsort(vals)

        s = input.xyz[idx[0]]
        d1 = input.xyz[idx[-1]]
        d2 = input.xyz[idx[-2]]

        # explicit output definitions
        sout = join(config['path'], 'output', 'nwchem', '%s_%s_%s_s.xyz' % (wildcards.id, wildcards.adduct, wildcards.cycle))
        d1out = join(config['path'], 'output', 'nwchem', '%s_%s_%s_d1.xyz' % (wildcards.id, wildcards.adduct, wildcards.cycle))
        d2out = join(config['path'], 'output', 'nwchem', '%s_%s_%s_d2.xyz' % (wildcards.id, wildcards.adduct, wildcards.cycle))

        shell('cp %s %s' % (s, sout))
        shell('cp %s %s' % (d1, d1out))
        shell('cp %s %s' % (d2, d2out))
