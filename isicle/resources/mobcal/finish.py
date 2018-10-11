import subprocess
from os.path import *
import pandas as pd
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW


def tail(f, lines=1):
    return [x.strip() for x in subprocess.check_output(['tail', '-n%s' % lines, f]).decode('ascii').split('\n')]


def boltzmann(infile, outfile):
    df = pd.read_csv(infile, sep='\t')

    g = df['DFT Energy'].values * 627.503
    mn = g.min()
    relG = g - mn
    b = np.exp(-relG / 0.5924847535)
    w = (b / b.sum()) * len(b)

    ws = DescrStatsW(df['Mean CCS'], weights=w, ddof=0)

    res = pd.Series([ws.mean, ws.std, ws.std_mean, ws.var, len(df.index)],
                    index=['mean', 'std', 'std_mean', 'var', 'N'])

    res.to_csv(outfile, sep='\t', index=False, header=True)


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
        return [m_mn, ccs_mn, ccs_std]
    else:
        return None
