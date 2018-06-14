# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 14:26:44 2015

@author: Dennis Thomas
"""
import glob
import logging
import re

from timeit import default_timer as timer

logger = logging.getLogger(__name__)


class AtomMassCharge:
    """
    This class contains the variables and methods to store, read and write 
    the mass and charges of all the atoms of a molecule. The file to read 
    from is the NWChem output file.

    Class is initialized with 4 arguments: molecule id (molecule_id),
    directory path to the NWChem output file (filedir), 
    the file name (filename), and the type of charges (chg_type;
    Lowdin vs. Mulliken)

    """

    def __init__(self, molecule_id, filedir, filename, chg_type='Lowdin_Charges'):
        self.molecule_id = molecule_id
        self.filedir = filedir
        self.filename = filename
        self.atom_mass = []
        self.atom_type = []
        self.atom_label = []
        self.atom_charge = []
        self.chg_type = chg_type

    # Function to read the atomic mass of each atom type from the
    # NWChem output file
    def read_atomMasses_fromNWChemOutputFile(self):
        """
        This function reads the atomic mass of each atom type in the 
        molecule. 

        INPUT: just the class member variables    
        OUTPUT: atom type, atomic mass
        USAGE:
           read_atomMasses_fromNWChemOutputFile()

        """
        line1 = 'Atomic Mass'
        line2 = 'Effective nuclear repulsion energy (a.u.)'

        fname = self.filedir + '/' + self.filename

        line_num1 = get_line_num(line1, fname)
        line_num2 = get_line_num(line2, fname)

        print line_num1
        print line_num2
        stop_read = 0

        fout = self.molecule_id + '_mass_charge.txt'
        g = open(fout, 'w')

        f = open(fname, 'r')
        line_counter = 0

        while stop_read == 0:
            line_counter += 1
            line = f.readline()
            if line_counter == line_num1:
                g.write(line)
            if (line_counter > line_num1 and
                    line_counter <= line_num2):
                g.write(line)
                line = line.strip()
                line_strs = line.split()
                if len(line_strs) == 2:
                    self.atom_type.append(line_strs[0])
                    self.atom_mass.append(float(line_strs[1]))
            if line_counter > line_num2:
                stop_read = 1

        f.close()
        g.write('\n')
        g.close()

    # Function to read the charges of each atom of the molecule from the
    # NWChem output file
    def read_atomCharges_fromNWChemOutputFile(self):
        """
        This function reads the charges of all the atoms in the 
        molecule. 

        INPUT: class member variables    
        OUTPUT: atom labels, atomic charge
        USAGE:
           read_atomCharges_fromNWChemOutputFile()

        """
        line1 = 'Total Density - Lowdin Population Analysis'
        line2 = 'Total Density - Mulliken Population Analysis'
        line3 = 'Atom       Charge   Shell Charges'

        search_line = ''
        if self.chg_type == 'Lowdin_Charges':
            search_line = line1
        elif self.chg_type == 'Mulliken_Charges':
            search_line = line2
        else:
            print 'Charge type cannot recognized.'

        print 'CHECKING IF THIS LINE WRITES'
        # Get the line number of the last line containing
        # the seaerch_line string, in the file.
        fname = self.filedir + '/' + self.filename
        f = open(fname, 'r')
        line_counter = 0
        search_line_num = 0
        stop_read = 0
        for line in f:
            line = line.strip()
            line_counter += 1
            if search_line in line:
                search_line_num = line_counter
        f.close()

        print search_line_num

        fout = self.molecule_id + '_mass_charge.txt'
        g = open(fout, 'a')

        f = open(fname, 'r')
        stop_read = 0
        line_counter = 0

        start_read = 0
        while stop_read == 0:
            line = f.readline()
            line_counter += 1
            if (line_counter >= search_line_num and
                    stop_read == 0):
                if line_counter == search_line_num:
                    g.write(line)
                rline = f.readline()
                line = rline.strip()
                if line3 in line:
                    g.write(rline)
                    start_read = 1
                    while start_read == 1:
                        line = f.readline()
                        g.write(line)
                        line = line.strip()

                        line_strs = line.split()
                        ncols = len(line_strs)
                        if ncols >= 4:
                            self.atom_label.append(line_strs[1])
                            self.atom_charge.append(float(line_strs[3]))
                        if not line:
                            start_read = 0
                            stop_read = 1
        g.close()
        f.close()

# Function to get the line number of a line containing
    # a string, in a file.


def get_line_num(linestr, infile):
    f = open(infile, 'r')
    print infile
    print linestr
    line_counter = 0
    for line in f:
        line = line.strip()  # remove '\n' at end of line
        line_counter += 1
        if linestr in line:
            f.close()
            return line_counter

#----------------------------------------------------------------------------


def createXYZFile_FromNWChemOutputFile(nwchem_dftoutputfile_dir, xyz_dir):
    """
        This function reads the xyz data from NWChem DFT output files, and 
        writes them out . 

        INPUT: NWChem DFT output file directory (nwchem_dftoutputfile_dir), 
                XYZ file directory (xyz_dir)
        OUTPUT: .xyz files in directory, "xyz_dir"
        USAGE:
           createXYZFile_FromNWChemOutputFile(<dft output dir>, < xyz dir>)

        """
  #  search_line1 = 'Optimization converged'  # use this if NWChem DFT output
 #   file were created during calculation.
    search_line1 = 'auto-z'
    search_line2 = 'Geometry'

    stop_line = 'Atomic Mass'
    flist = nwchem_dftoutputfile_dir + '/*'
    print flist

    for mfile in glob.glob(flist):

        x = []
        y = []
        z = []
        atom = []
        numatoms = 0
        f = open(mfile, 'r')
        line_counter = 0
        search_line1_num = 0
        start_read = 0
        find_search_line2 = 0
        start_skipping = 0
        stop_read = 0
        for line in f:
            line = line.strip()
            line_counter += 1
            if search_line1 in line:
                search_line1_num = line_counter
                find_search_line2 = 1
            if (find_search_line2 == 1 and search_line2 in line):

                start_skipping = 1
            if(start_skipping == 1):
                start_read = start_read + 1

            if (start_read >= 8 and stop_read == 0):

                strs = line.split()
                if (len(strs) > 2):
                    atom.append(strs[1])
                    x.append(strs[3])
                    y.append(strs[4])
                    z.append(strs[5])
                    numatoms = numatoms + 1
                   # print atom[numatoms-1], x[numatoms-1], y[numatoms-1], z[numatoms-1]
                if stop_line in line:
                    stop_read = 1
                    # print numatoms

        f.close()

        # write xyz file
        s1 = re.split("/|,|\\\|\n", mfile)

        strs = s1[len(s1) - 1].split('.')
        fname = xyz_dir + "/" + strs[0] + '.xyz'

        fout = open(fname, 'w')
        fout.write(str(numatoms) + '\n')
    #    fout.write('Geometry' + '\n')
        fout.write('\n')

        for a, xi, yi, zi in zip(atom, x, y, z):
            fout.write("{}\t{}\t{}\t{}\n".format(a, xi, yi, zi))
        fout.close()

        msg = "[nwchem]: Read xyz values from file " + mfile
        logger.info(msg)
        msg = "[nwchem]:search line 1 no. = " + str(search_line1_num)
        logger.info(msg)
