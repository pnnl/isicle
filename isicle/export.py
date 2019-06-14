import glob
from os.path import *
import pandas as pd
import os
import shutil


def gather_inputs():
    # input files
    files_in = glob.glob('input/*.smi')
    files_in.extend(glob.glob('input/*.inchi'))
    if len(files_in) < 1:
        raise IndexError('No input files found')

    # parse input file content
    d = {}
    d['Index'] = [splitext(basename(x))[0] for x in files_in]
    d['Format'] = [splitext(basename(x))[1][1:] for x in files_in]
    d['String'] = []

    for file_in in files_in:
        with open(file_in, 'r') as f:
            d['String'].append(f.readlines()[0].strip())

    # input dataframe
    return pd.DataFrame(d)


def ccs(path, mode='standard'):
    # inputs
    df_in = gather_inputs()

    # standard
    if mode == 'standard':
        # output files
        files_out = glob.glob('output/mobility/mobcal/calibrated_ccs/*.tsv')
        if len(files_out) < 1:
            raise IndexError('No output files found')

        # parse output file content
        d = {}
        d['Index'] = [splitext(basename(x))[0].split('_')[0] for x in files_out]
        d['Adduct'] = [splitext(basename(x))[0].split('_')[1] for x in files_out]

        content_out = pd.concat([pd.read_csv(x, sep='\t') for x in files_out], axis=0, ignore_index=True)
        d['CCS'] = content_out['ccs'].values
        d['CCS_std'] = content_out['ccs_std'].values
        d['N'] = content_out['n'].values

        # output dataframe
        df_out = pd.DataFrame(d)

        # join
        df = df_in.merge(df_out, on=['Index'], how='left')

        # save
        df.to_csv(path, sep='\t', index=False)

    # lite
    elif mode == 'lite':
        # output files
        files_out = glob.glob('output/mobility/impact/ccs/*.ccs')
        if len(files_out) < 1:
            raise IndexError('No output files found')

        # parse output file content
        d = {}
        d['Index'] = [splitext(splitext(basename(x))[0])[0].split('_')[0] for x in files_out]
        d['Adduct'] = [splitext(splitext(basename(x))[0])[0].split('_')[1] for x in files_out]
        d['Buffer Gas'] = [splitext(splitext(basename(x))[0])[1][1:] for x in files_out]
        d['CCS'] = []

        for file_out in files_out:
            with open(file_out, 'r') as f:
                d['CCS'].append(f.readlines()[0].strip())

        # output dataframe
        df_out = pd.DataFrame(d)

        # join
        df = df_in.merge(df_out, on=['Index'], how='left')

        # save
        df.to_csv(path, sep='\t', index=False)

    # other
    else:
        raise ValueError('Mode must be "lite" or "standard"')


def shifts(path):
    # inputs
    df_in = gather_inputs()

    # output files
    files_out = glob.glob('output/shifts/*.tsv')
    if len(files_out) < 1:
        raise IndexError('No output files found')

    # output directory
    if not exists(path):
        os.makedirs(path)

    # save index
    df_in.to_csv(join(path, 'index.tsv'), sep='\t', index=False)

    # copy results
    [shutil.copy2(x, path) for x in files_out]


def energy(path):
    # output files
    files_out = glob.glob('output/energy/*.energy')
    if len(files_out) < 1:
        raise IndexError('No output files found')

    # parse output file content
    d = {}
    d['Index'] = [splitext(splitext(basename(x))[0])[0].split('_')[0] for x in files_out]
    d['Total DFT Energy'] = []

    for file_out in files_out:
        with open(file_out, 'r') as f:
            d['Total DFT Energy'].append(f.readlines()[0].strip())

    # output dataframe
    df_out = pd.DataFrame(d)

    # join
    df = df_in.merge(df_out, on=['Index'], how='left')

    # save
    df.to_csv(path, sep='\t', index=False)
