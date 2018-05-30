# -*- coding: utf-8 -*-
"""
Processes IMPACT files

@author: MontyPypi
"""

import os

path = r'C:\Users\nune558\Documents\CCS Calc\DSSTox'
for adduct in ['deprotonated', 'protonated', 'sodiated']:

    # Locate and open IMPACT file
    fname = r'%s\IMPACT' % adduct
    full_fname = path + fname + '.txt'
    if os.path.exists(full_fname):
        f = open(full_fname)
        flines = f.readlines()[8:]
        f.close()
        
        # Process and save in processed.txt file
        f = open(path + fname + '_processed.txt', 'w')
        for line in flines:
            tokens = [token for token in line.split(' ') if len(token) > 0]
            molid = line.split('\\')[-1].split('+')[0]
            ccs = tokens[-2]
            f.write('\t'.join([molid, ccs]))
            f.write('\n')
        f.close()
