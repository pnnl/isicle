# -*- coding: utf-8 -*-\
"""
Created on Sat Oct 24 23:23:44 2015

@author: Dennis Thomas
"""
import nwchem
import geometry

# Function to create the the input files for MOBCAl (.mfj and mobcal.in file)
class MobcalInputFiles:
    """
    This class contains the variables and methods to create the MOBCAL input
    files (.mfj and mobcal.in).
    
    """
    def __init__(self, molecule_id, mfj_filedir, mfj_filename, nconf = 1):
        self.mfj_filedir = mfj_filedir
        self.mfj_filename = mfj_filename
        self.molecule_id = molecule_id
        self.xyz_filedir = ''
        self.xyz_filename = ''
        self.nwc_filedir = ''
        self.nwc_filename = ''
        self.chg_type = ''
        self.natoms = 0
        self.nconformers = nconf # number of conformers
        self.data = []
        
    def create_mfjFile_fromXYZ_NWChemOutput(self, xyz_filedir, xyz_filename,
                 nwc_filedir, nwc_filename, chg_type = 'Lowdin_Charges'):
        
        self.xyz_filedir = xyz_filedir
        self.xyz_filename = xyz_filename
        self.nwc_filedir = nwc_filedir
        self.nwc_filename = nwc_filename
        self.chg_type = chg_type
        
        # Read the geometry (x,y,z) atom coordinates of the molecule
        d = geometry.XYZData(self.xyz_filedir,self.xyz_filename)
        d.read_xyz()
        # Write out the data (just for checking)
        d.write_xyz('./','out_'+self.xyz_filename)

        # Get the mass and charges of each atom in the molecule
        c = nwchem.AtomMassCharge(self.molecule_id, self.nwc_filedir,
                                            self.nwc_filename,'Lowdin_Charges')
        c.read_atomMasses_fromNWChemOutputFile()
        c.read_atomCharges_fromNWChemOutputFile()
        
        self.natoms = len(c.atom_label)
        #Write the .mfj file
        fname = self.mfj_filedir + './' + self.mfj_filename
        f = open(fname,'w')
        
        f.write(str(self.nconformers) + '\n')
        f.write(str(self.natoms) + '\n')
        f.write('ang\n')
        f.write('calc\n')
        f.write('1.000\n')
  
        masses = []
        for i in range(self.natoms):
           # x = d.x[i]
           # y = d.y[i]
           # z = d.z[i]
            #chg = c.atom_charge[i]
            stop = 0
            j = 0
            while stop == 0:
                if c.atom_type[j] == c.atom_label[i]:
                    v1 = c.atom_mass[j]
                   
                    m = round(v1,1) 
                    masses.append(int(m))
                    stop = 1
                else:
                    j += 1
        
        for x, y, z, m, chg in zip(d.x, d.y, d.z, masses, c.atom_charge):
            f.write("{}\t\t{}\t\t{}\t\t{}\t\t{}\n".format(str(x), str(y), 
                    str(z), str(m), str(chg)))
            
        f.close()
    
    # Function to create mobcal.in file
    def create_mobcalin_file(self):
        fname = 'mobcal.in'
        f = open(fname,'w')
        f.write(self.mfj_filename + '\n')
        f.write(self.molecule_id + '.out' + '\n')
        f.write('5013489\n')
        f.close()