from os.path import *
from resources.utils import *
import shutil
import subprocess
import os
from string import Template
import glob

# snakemake configuration
configfile: 'config.yaml'

IDS, ADDUCTS = glob_wildcards(join(config['path'], 'input', '{id}_{adduct}.mol2'))


def cycles():
    return ['%03d' % x for x in range(config['cycles'] + 1)]


def frames():
    return ['%03d' % x for x in range(config['nframes'])]


# a pseudo-rule that collects the target files
rule all:
    input:
        expand(join(config['path'], 'output', 'selected', 'xyz', '{id}_{adduct}_{cycle}_s.xyz'),
               id=IDS, adduct=ADDUCTS, cycle=cycles()[1:])

rule antechamber:
    input:
        mol2 = join(config['path'], 'input', '{id}_{adduct}.mol2')
    output:
        tmp = join(config['path'], 'output', 'antechamber', 'tmp', '{id}_{adduct}', '{id}_{adduct}.input.mol2'),
        ac = join(config['path'], 'output', 'antechamber', 'tmp', '{id}_{adduct}', '{id}_{adduct}.output.mol2'),
        frcmod = join(config['path'], 'output', 'antechamber', 'frcmod', '{id}_{adduct}.frcmod'),
        mol2 = join(config['path'], 'output', 'antechamber', 'mol2', '{id}_{adduct}.mol2')
    log:
        ac = join(config['path'], 'output', 'antechamber', 'logs', '{id}_{adduct}.antechamber.log'),
        parmchk = join(config['path'], 'output', 'antechamber', 'logs', '{id}_{adduct}.parmchk2.log')
    run:
        mol = read_mol(input.mol2)

        # this doesn't always work
        ###################################
        charge = mol.total_partial_charge() - 1
        ###################################

        natoms = mol.natoms()

        if wildcards.adduct == '+Na':
            idx, content = pop_atom(input.mol2, output.tmp, atom='Na')
            charge -= len(idx)
        else:
            shutil.copy2(input.mol2, output.tmp)

        # abspath
        tmp = abspath(output.tmp)
        ac = abspath(output.ac)
        frcmod = abspath(output.frcmod)
        log_ac = abspath(log.ac)
        log_parmchk = abspath(log.parmchk)

        cwd = os.getcwd()
        os.chdir(join(config['path'], 'output', 'antechamber', 'tmp', '%s_%s' % (wildcards.id, wildcards.adduct)))

        subprocess.call('antechamber -i %s -fi mol2 -o %s -fo mol2 -c bcc -s -du -nc %.4f &> %s' %
                        (tmp, ac, charge, log_ac), shell=True)
        subprocess.call('parmchk2 -i ANTECHAMBER_AC.AC -f ac -o %s &> %s' %
                        (frcmod, log_parmchk), shell=True)
        os.chdir(cwd)

        if wildcards.adduct == '+Na':
            push_atom(output.ac, output.mol2, idx, content)
        else:
            shutil.copy2(output.ac, output.mol2)

rule tleap:
    input:
        mol2 = rules.antechamber.output.mol2,
        frcmod = rules.antechamber.output.frcmod
    output:
        config = join(config['path'], 'output', 'tleap', 'config', '{id}_{adduct}.config'),
        prmtop = join(config['path'], 'output', 'tleap', 'prmtop', '{id}_{adduct}.top'),
        inpcrd = join(config['path'], 'output', 'tleap', 'inpcrd', '{id}_{adduct}.crd')
    log:
        join(config['path'], 'output', 'tleap', 'logs', '{id}_{adduct}.log')
    run:
        with open('resources/amber/tleap.template', 'r') as f:
            t = Template(f.read())

        d = {'frcmod': input.frcmod,
             'mol2': input.mol2,
             'prmtop': join(config['path'], 'output', 'tleap', 'prmtop', '%s_%s.top' % (wildcards.id, wildcards.adduct)),
             'inpcrd': join(config['path'], 'output', 'tleap', 'inpcrd', '%s_%s.crd' % (wildcards.id, wildcards.adduct))}

        with open(output.config, 'w') as f:
            f.write(t.substitute(d))

        shell('tleap -s -f {output.config} &> {log}')

rule sander_em:
    input:
        mol2 = rules.antechamber.output.mol2,
        prmtop = rules.tleap.output.prmtop,
        inpcrd = rules.tleap.output.inpcrd
    output:
        config = join(config['path'], 'output', 'sander', 'em', '{id}_{adduct}.mdin'),
        rst = join(config['path'], 'output', 'sander', 'em', '{id}_{adduct}.rst'),
        out = join(config['path'], 'output', 'sander', 'em', '{id}_{adduct}.out')
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

        shell('sander -O -i {output.config} -o {output.out} -c {input.inpcrd} -p {input.prmtop} -r {output.rst}')

rule sander:
    input:
        mol2 = rules.antechamber.output.mol2,
        rst = rules.sander_em.output.rst,
        prmtop = rules.tleap.output.prmtop
    output:
        config = expand(join(config['path'], 'output', 'sander', 'anneal', '{cycle}', '{{id}}_{{adduct}}.mdin'), cycle=cycles()),
        rst = expand(join(config['path'], 'output', 'sander', 'anneal', '{cycle}', '{{id}}_{{adduct}}.rst'), cycle=cycles()),
        crd = expand(join(config['path'], 'output', 'sander', 'anneal', '{cycle}', '{{id}}_{{adduct}}.crd'), cycle=cycles()),
        out = expand(join(config['path'], 'output', 'sander', 'anneal', '{cycle}', '{{id}}_{{adduct}}.out'), cycle=cycles())
    run:
        # iterate SA steps
        for i in range(config['cycles'] + 1):

            # explicit output definitions
            cfg = join(config['path'], 'output', 'sander', 'anneal', '%03d' % i, '%s_%s.mdin' % (wildcards.id, wildcards.adduct))
            rst_out = join(config['path'], 'output', 'sander', 'anneal', '%03d' % i, '%s_%s.rst' % (wildcards.id, wildcards.adduct))
            crd = join(config['path'], 'output', 'sander', 'anneal', '%03d' % i, '%s_%s.crd' % (wildcards.id, wildcards.adduct))
            out = join(config['path'], 'output', 'sander', 'anneal', '%03d' % i, '%s_%s.out' % (wildcards.id, wildcards.adduct))

            # first SA step
            if i == 0:
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

                with open(cfg, 'w') as f:
                    f.write(t.substitute(d))

                rst_in = input.rst

            # remaining SA steps
            else:
                with open('resources/amber/sander_anneal.template', 'r') as f:
                    t = Template(f.read())

                d = {'mol2': input.mol2}

                with open(cfg, 'w') as f:
                    f.write(t.substitute(d))

                rst_in = join(config['path'], 'output', 'sander', 'anneal', '%03d' % (i - 1), '%s_%s.rst' % (wildcards.id, wildcards.adduct))

            cmd = 'sander -O -p %s -c %s -i %s -o %s -r %s -x %s' % (input.prmtop, rst_in, cfg, out, rst_out, crd)
            shell(cmd)

rule extract_frames:
    input:
        prmtop = rules.tleap.output.prmtop,
        out = join(config['path'], 'output', 'sander', 'anneal', '{cycle}', '{id}_{adduct}.out'),
        crd = join(config['path'], 'output', 'sander', 'anneal', '{cycle}', '{id}_{adduct}.crd')
    output:
        trajin = join(config['path'], 'output', 'sander', 'extracted', 'trajin', '{id}_{adduct}_{cycle}_{frame}.trajin'),
        mol2 = join(config['path'], 'output', 'sander', 'extracted', 'mol2', '{id}_{adduct}_{cycle}_{frame}.mol2'),
    log:
        join(config['path'], 'output', 'sander', 'extracted', 'logs', '{id}_{adduct}_{cycle}_{frame}.log')
    run:
        frame = select_frames(input.out,
                              nsamples=config['nframes'],
                              low=config['low'],
                              high=config['high'])[int(wildcards.frame)]

        shell('echo "trajin {input.crd} %s %s" > {output.trajin}' % (frame, frame))
        shell('echo "trajout {output.mol2} mol2" >> {output.trajin}')
        shell('cpptraj {input.prmtop} {output.trajin} > {log}')

rule convert:
    input:
        mol2a = rules.extract_frames.output.mol2,
        mol2b = rules.antechamber.input.mol2
    output:
        xyz = join(config['path'], 'output', 'sander', 'extracted', 'xyz', '{id}_{adduct}_{cycle}_{frame}.xyz')
    run:
        standardizeMol2(input.mol2a, input.mol2b, output.xyz)

rule calculate_rmsd:
    input:
        xyz = rules.convert.output.xyz
    output:
        rmsd = join(config['path'], 'output', 'selected', 'rmsd', '{id}_{adduct}_{cycle}_{frame}.rmsd')
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
        xyz = expand(join(config['path'], 'output', 'sander', 'extracted', 'xyz', '{{id}}_{{adduct}}_{{cycle}}_{frame}.xyz'), frame=frames()),
        rmsd = expand(join(config['path'], 'output', 'selected', 'rmsd', '{{id}}_{{adduct}}_{{cycle}}_{frame}.rmsd'), frame=frames())
    output:
        s = join(config['path'], 'output', 'selected', 'xyz', '{id}_{adduct}_{cycle}_s.xyz'),
        d1 = join(config['path'], 'output', 'selected', 'xyz', '{id}_{adduct}_{cycle}_d1.xyz'),
        d2 = join(config['path'], 'output', 'selected', 'xyz', '{id}_{adduct}_{cycle}_d2.xyz')
    run:
        vals = []
        for f in input.rmsd:
            vals.append(float(read_string(f)))

        vals = np.array(vals)
        idx = np.argsort(vals)

        s = input.xyz[idx[0]]
        d1 = input.xyz[idx[-1]]
        d2 = input.xyz[idx[-2]]

        shell('cp %s {output.s}' % s)
        shell('cp %s {output.d1}' % d1)
        shell('cp %s {output.d2}' % d2)
