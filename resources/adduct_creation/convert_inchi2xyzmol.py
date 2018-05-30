# -*- coding: utf-8 -*-
"""
Created on Sat Sep 03 01:00:25 2016

@author: Dennis Thomas
"""

import sys
sys.path.insert(0, '../')
sys.path.insert(0, '../../')

import geometry
import pka
import readio
import csv
import pybel
import openbabel
import datetime
import logging
import glob
import numpy


from timeit import default_timer as timer


def main():
    # read inchi file list

    # read input parameters from file, parameters.in

    dt = datetime.datetime.now()
    dtstr = dt.strftime("%A, %B %d, %Y %I:%M %p")

    # create a name for the log file based on start date and time

    logfile = 'run_' + str(dt.month) + "_" + str(dt.day) + "_" + str(dt.year) + "_"
    logfile = logfile + str(dt.hour) + "_" + str(dt.minute) + "_" + str(dt.second) + ".log"
    print logfile

    logging.basicConfig(filename=logfile, level=logging.INFO)
    logger = logging.getLogger(__name__)
    msg = "Started " + dtstr
    logger.info(msg)

    fname = 'parameters.in'
    msg = "Reading parameters from " + fname
    logger.info(msg)

    uio = readio.IOFile(fname)
    uio.readFile()

    # read inCHi list
    msg = "Reading InChI list"
    logger.info(msg)
    d = geometry.InChI(uio.inchilist_dir, uio.inchilist_file)
    d.read_inChiFile()

    # Convert inChI to .mol format
    time_inchi2molxyz_seconds = []
    mol3Dfiles = []
    for i in range(0, d.num_cmpnds):
        tstart = timer()

        inchi_str = d.inchi[i]
        cpd_id = d.object_id[i]

        (mol, molfile, mol2Dfile, mol3Dfile) = geometry.inchi2mol(inchi_str, cpd_id,
                                                                  "/", uio.inchi2mol_dir, uio.forcefield)

    # Create .xyz. file
#        geometry.createXYZFile(mol,mol3Dfile,uio.inchi2xyz_dir)
        tend = timer()
        mol3Dfiles.append(mol3Dfile)
        time_inchi2molxyz_seconds.append(tend - tstart)

    with open('time_inch2molxyz.csv', 'wb') as fp:
        a = csv.writer(fp, delimiter=',')
        a.writerow(["mol_id", "mol3D_file", "time_seconds"])

        for i in range(0, d.num_cmpnds):
            out = [d.object_id[i], mol3Dfiles[i], str(time_inchi2molxyz_seconds[i])]
            a.writerow(out)

    fp.close()


if __name__ == '__main__':
    main()
