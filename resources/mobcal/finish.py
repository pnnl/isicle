import subprocess
from os.path import *
import pandas as pd
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW


def tail(f, lines=1):
    return [x.strip() for x in subprocess.check_output(['tail', '-n%s' % lines, f]).decode('ascii').split('\n')]


def _func(grp):
    ref = grp['G'].min()
    grp['RelG'] = grp['G'] - ref
    grp['B'] = np.exp(-grp['RelG'] / 0.5924847535)
    grp['Weights'] = (grp['B'] / grp['B'].sum()) * len(grp.index)

    ws = DescrStatsW(grp['Mean CCS'], weights=grp['Weights'], ddof=0)

    return pd.Series([ws.mean, ws.std, ws.std_mean, ws.var, len(grp.index)],
                     index=['mean', 'std', 'std_mean', 'var', 'N'])


def boltzmann(infile, outfile):
    df = pd.read_csv(infile, sep='\t')

    # df['File'] = [x[-1] for x in df['File'].str.rsplit('/', 1).tolist()]

    info = [x[1].split('+') for x in df['File'].str.split('_').tolist()]

    df['ID'] = ['molid' + str(x[0]) for x in info]
    df['Adduct'] = ['+' + ''.join(i for i in x[1] if not i.isdigit()) for x in info]

    df['G'] = df['DFT Energy'] * 627.503
    df = df.groupby(by=['ID', 'Adduct']).apply(_func).reset_index()

    df.to_csv(outfile, sep='\t', index=False)


def parse_mobcal(f):
    done = False

    lines = tail(f, lines=10)

    for line in lines:
        if "average (second order) TM mobility" in line:
            m_mn = float(line.split('=')[-1])
        elif "average TM cross section" in line:
            ccs_mn = float(line.split('=')[-1])
        elif "standard deviation" in line:
            ccs_std = float(line.split('=')[-1])
            done = True
    if done is True:
        return [basename(f), m_mn, ccs_mn, ccs_std]
    else:
        return None
