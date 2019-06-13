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
    ikey_in = [splitext(basename(x))[0] for x in files_in]
    fmt_in = [splitext(basename(x))[1][1:] for x in files_in]
    content_in = []
    for file_in in files_in:
        with open(file_in, 'r') as f:
            content_in.append(f.readlines()[0].strip())

    # input dataframe
    return pd.DataFrame({'InChI Key': ikey_in, 'String': content_in, 'Format': fmt_in})


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
        ikey_out = [splitext(basename(x))[0].split('_')[0] for x in files_out]
        adduct_out = [splitext(basename(x))[0].split('_')[1] for x in files_out]
        content_out = pd.concat([pd.read_csv(x, sep='\t') for x in files_out], axis=0, ignore_index=True)

        # output dataframe
        df_out = pd.DataFrame({'InChI Key': ikey_out,
                               'Adduct': adduct_out,
                               'CCS': content_out['ccs'].values,
                               'CCS_std': content_out['ccs_std'].values,
                               'N': content_out['n'].values})

        # join
        df = df_in.merge(df_out, on=['InChI Key'], how='left')

        # save
        df.to_csv(path, sep='\t', index=False)

    # lite
    elif mode == 'lite':
        raise ValueError('Lite not yet supported')

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
    raise ValueError('Free energy module not yet supported')
