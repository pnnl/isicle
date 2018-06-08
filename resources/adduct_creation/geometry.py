import csv
import pybel
import sys
import logging
import openbabel
import numpy
import glob


logger = logging.getLogger(__name__)

class Atom:
    def __init__(self,x,y,z,r,chg,atomt):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.chg = chg
        self.atomtype = atomt
    
    
class XYZData:
    """
    This class contains the variables and methods to store, read and write 
    (x,y,z) coordinates of atoms in a molecule. Files to read adn write are in
    .xyz format.
    
    """
    def __init__(self,filedir,filename):
        self.filedir = filedir
        self.filename  = filename
        self.molecule_id = filename.split('.')[0]
        self.atom = []
        self.x = []
        self.y = []
        self.z = []
        self.numatoms = []
        
     # Function to read in the x,y,z data   
    def read_xyz(self):
        """
        This function reads the (x,y,z) coordinates of each atom and 
        its type (e.g., 'O' for oxygen atom, 'H' for hydrogen atom, etc.).
        
        INPUT: path to file directory(filedir) and file name (filename)    
        OUTPUT: atom type, (x,y,z) coordinates
        USAGE:
            read_xyz(<path-to-file-directory>,<file-name>)
        EXAMPLE:
            read_xyz('./xyz_files/','test.xyz')
        XYZ FILE FORMAT:
            four-column data in the order: 
                <atom-type> <x-coord> <y-coord> <z-coord>
            Coordinates unit: in Angstrom
        """
        print "The file - " + self.filename + " - is located in directory:"
        print self.filedir
     
        fname = self.filedir + "/" + self.filename
        
        print fname
        f = open(fname, 'r')
        
        self.numatoms = 0
        for line in f:
            values = line.split()
            
            if len(values) == 4:
                self.numatoms += 1
                values = line.split()
                self.atom.append(values[0])
                self.x.append(values[1])
                self.y.append(values[2])
                self.z.append(values[3])
        f.close()
        
    # Function to write out .xyz files
    def write_xyz(self,filedir,filename):
        """
        This function writes out the (x,y,z) coordinates of each atom and 
        its type (e.g., 'O' for oxygen atom, 'H' for hydrogen atom, etc.), in 
        .xyz file format.
        
        INPUT: path to file directory(filedir) and file name (filename)    
        OUTPUT: atom type, (x,y,z) coordinates
        USAGE:
            write_xyz(<path-to-file-directory>,<file-name>)
        EXAMPLE:
            write_xyz('./','test_out.xyz')
        XYZ FILE FORMAT:
            four-column data in the order: 
                <atom-type> <x-coord> <y-coord> <z-coord>
            Coordinates unit: in Angstrom
        """
        fname = filedir + "/" + filename
        
        f = open(fname, 'w')
        f.write(str(self.numatoms)+ '\n')
        f.write('Geometry' + '\n')
        f.write('\n')
        
        for atom, x, y, z in zip(self.atom, self.x, self.y, self.z):
            f.write("{}\t{}\t{}\t{}\n".format(atom,str(x),str(y),str(z)))
        f.close()

    # Function to create XYZ grid

#==============================================================================
#     # Function to create new adducts from the standard molecule
#     def addAtom(self,nearatom_num,newatom_radius,newatom_type):
#         """
#             This function:
#             1) adds adds a new atom 
#             2) saves the new molecule in an .xyz file.
#             
#             INPUT:
#                 adduct_type:  '+H', '-H','+Na'
#                 
#         """
#         nearatom_type = self.atom(nearatom_num - 1)
#         
#         if newatom_type == 'H' or newatom_type == 'Na':
#             xi = self.x(nearatom_num - 1)
#             yi = self.y(nearatom_num - 1)
#             zi = self.z(nearatom_num - 1)
#             
#             nearatom = Atom(xi,yi,zi,r,chg,nearatom_type)
#             [xj,yj,zj] = findPosition()
#             
#             #new atom 
#             newatom = openbabel.OBAtom
#             newatom.OBAtom.SetVector(xj,yj,zj)
#             
#             
#         else:
#             msg = "[geometry.py(XYZData.addAtom)]: Cannot add %s."
#             msg = "Only hydrogen or sodium can be added."
#             msg % (atom_type)
#             sys.exit(msg)  
#         
#     
#         
#         # remove Hydrogen
#     def removeAtom(self,atom_num):
#         """
#             This function:
#             1) removes hydrogen atom, given its ranking (atom_num) in the
#                .xyz file. Updates the number of atoms.
#             
#             INPUT:
#                 atom_num:  rank of the hydrogen atom, as listed in the .xyz
#                            file.
#             
#         
#         """
#         if self.atom[atom_num -1] == 'H':
#             
#             del self.x[atom_num - 1]
#             del self.y[atom_num - 1]
#             del self.z[atom_num - 1]
#             del self.atom[atom_num - 1]
#             self.numatoms = self.numatoms - 1
#             msg = "[geometry.py(XYZData.removeAtom)]: removing %s"
#             msg = msg + ", located at %d, %d, and %d."
#             msg % (self.atom[atom_num-1],self.x[atom_num - 1],
#             self.y[atom_num - 1], self.z[atom_num - 1])
#             
#             logging.info(msg)
#         else:
#             msg = "[geometry.py(XYZData.removeAtom)]: The atom to remove"  
#             msg = msg + " from the molecule is not hydrogen."
#             sys.exit(msg)
#     
#==============================================================================
        
            

#read xyz file of the molecule
# select atom line number and type 
# to remove H, we need to know which one it is.
# to add H, we need to know to which atom we want to associate H with.
# to add Na, we need to know  to which atom we want to associate Na with.

    
# class InChI
class InChI:
    """
    This class reads and stores the InChI strings of compounds from a CSV file.
    """
    def __init__(self,filedir,filename):
        self.filedir = filedir
        self.filename = filename
        self.num_cmpnds = []
        self.object_id = [] # compound id
        self.cmpnd = [] # compound name
        self.inchi = [] # InChI string
        self.inchi_key = [] # InChI key
        self.mw = [] #  molecular weight
        self.adduct_type = [] # adduct type
        self.src_molid = [] # source mol id
       
    
    # Function to read the InChI from a tab-delimited .txt file or .csv file.
    def read_inChiFile(self):
        """
        This function reads InChI strings from a tab-delimited .txt file. 
        
        INPUT: path to file directory(filedir) and file name (filename)    
        OUTPUT: 
        USAGE:
            read_inChiFile()
        EXAMPLE:
            read_inChiFile()
        FILE COLUMN FORMAT (tab-delimited .txt): 5 columns
          Object ID, Compound,InChI,InChI-Key, Monoisotopic-Molecular-Weight
          
        """
        print "The file - " + self.filename + " - is located in directory:"
        print self.filedir
     
        fname = self.filedir + "/" + self.filename

        if '.csv' in self.filename:
           
            f = csv.reader(open(fname),delimiter = ',')
        elif '.txt' in self.filename:
            f = csv.reader(open(fname),delimiter = '\t')
#        if '.txt' not in self.filename:
#        
#           err_msg = ("Error: The InChi list has to be in a 5-column " + 
#           "tab-delimted .txt file.") 
#           sys.exit(err_msg)

        print fname

      #  f = open(fname, 'r')
        line_num = 0
        for values in f:
                  
            line_num += 1
            if(line_num ==1):
                index0 = values.index('Object ID') 
                index1 = values.index('Compound')
                index2 = values.index('InChI')
                index3 = values.index('InChI-Key')
                index4 = values.index('Monoisotopic-Molecular-Weight')
                index5 = values.index('Adduct Type')
                index6 =  values.index('Source_MolID')
                    
            if line_num > 1:
               
               
                self.object_id.append(values[index0])
                self.cmpnd.append(values[index1])
                self.inchi.append(values[index2])
                self.inchi_key.append(values[index3])
                self.mw.append(values[index4])
                self.adduct_type.append(values[index5])
      
        self.num_cmpnds = line_num  - 1
	
       # f.close()      
    
# Function to convert InChI string to mol format using the pybel module

def inchi2mol(inchi_str,cpd_id,file_prefix,filedir,ffield):
    """
    This function:
    1) Creates a molecule from an inChI string using the 
    pybel.readstring function.
    2) Creates a 2D representation of
    the molecule, and saves it as an image file (in .png format)
    3) Creates a 3D coordinate file in .mol format.
    INPUT:
        inchi_str:  inChI string
        inchi_key: inChI key
        file_prefix:  The prefix for the .mol file
        filedir: The directory where the .mol file is saved.
    
    OUTPUT:
        mol: The molecule created from the inChI string
        molfile: The name of the .mol file.
        mol2Dfile:  The name of the image file (in .png format)
        
    Note: The files are saved in the format {mol}<id_><inChI_key><file ext>
    """
    #strs = inchi_key.split("InChIKey=")

    print inchi_str
    try:
        mol = pybel.readstring("inchi",inchi_str)
        mol.addh() # not necessary, because pybel make3D will add hydrogen
        
        molfile =filedir + file_prefix + cpd_id + ".mol"
        
        mol.write('mol',molfile,True)
        mol2Dfile = filedir + file_prefix + cpd_id + "_2D.png"
        mol.draw(show=False,filename=mol2Dfile,usecoords=False,update=False)
        print "molecule saved in .mol file format as " + molfile
#        print "2D representation of the molecule saved as " + mol2Dfile
#        print mol2dfile
     
        # Optimize 3D geometry of the molecule using pybel module's make3D() 
        # function    
        mol3Dfile = filedir + file_prefix + cpd_id + "_3D.mol"
        mol.make3D(forcefield=ffield,steps=50)
        mol.localopt(forcefield=ffield,steps=500)
       # mol.draw(show=False,filename=mol3Dfile,usecoords=True,update=True)
    
#        print "Created and optimized 3D geometry of the molecule"
#        print "forcefield."
        
            
        mol.write('mol',mol3Dfile,True)
        del mol
        
        return (None, molfile,mol2Dfile,mol3Dfile)
    except IOError:
        return (None, None, None, None)

    
# Function to create xyz file from optimized 3D geometry of the molecule
def createXYZFile(mol,mol3Dfile,filedir):
    """
    This function creates a .xyz file from the molecule 3D coordinate file in mol format.
    INPUT:
        mol:  pybel molecule object 
        mol3Dfile:  .mol file
        filedir: The directory where the .xyz file is saved.
    
    OUTPUT:
        xyzfile: The .xyz file
    """
    print "mol3Dfile = " + mol3Dfile
    strs0 = mol3Dfile.split("_3D.mol")
    strs1 = strs0[0].split("/")
    xyzfile = filedir + '/' + strs1[-1] + ".xyz"
    print "xyz file name: " + xyzfile
    
    mol.write("xyz",xyzfile,True)
    

#Function to add atom to molecule in .mol format
def addAtomToMol(mol,add_element,near_atomnum,covalent_bond_type,bonded):
    
    # get the van der Waals radius of the atoms in the current molecule
    # and create the x,y,z,radius array
    """
    This function adds an atom (add_element) to a molecule (mol) at the site
    near an atom (near_atomnum). It uses the covalent radius of the atoms
    to determine the accessible site for placing the the new atom.
    
    INPUT:
        mol -           pybel molecule object
        add_element -   element symbol (string) of the atom to be added
        near_atomnum -  atom number of the atom near which the new atom is 
                        placed
        bond_type:  bond order (current value used is "single" (only for covalent bond; otherwise, it's not added)
        bonded - ("yes" for covalent; 'no' for non-bonded (van der Waal))
    OUTPUT:
        mol - pybel molecule object containing the new atom
        xyzr - molecule with atom (x,y,z) coordinates and radius
    
    USAGE:
        (<new pybel mol object>,<xyzr matrix>) = 
        addAtomToMol(<pybel mol object>,<element symbol>,<atom number>)
        
        e.g., (molout,xyzr) = addAtomToMol(mol,'H','3')
    """
    
    etab = openbabel.OBElementTable()     # element object
     
    atoms_radius = []
    xyzr = []

    for atom in mol.atoms:
        atomic_num = atom.atomicnum # get atomic number of each atom
       # vdw_rad = etab.GetVdwRad(atomic_num)
        if(bonded == "yes"): # if bond type is covalent bond
            # use covalent radius            
            atom_rad = etab.GetCovalentRad(atomic_num)
            

        elif(bonded == "no"): # bond type is non-covalent bond
            # use van der Waals radius (but for now use covalent)
           #atom_rad = etab.GetCovalentRad(atomic_num) 
           atom_rad = etab.GetVdwRad(atomic_num) 
        atoms_radius.append(atom_rad)
        
        a = atom.coords + (atom_rad,)  # create x,y,z,radius array
        xyzr.append(a)
       
    # get the atomic number and covalent radius of the new atom to be added 
    
    newatom_atomic_num = etab.GetAtomicNum(add_element) # get atomic number
    if(bonded == "yes"): # if bond type is covalent bond
        #newatom_rad = etab.GetCovalentRad(newatom_atomic_num) # get covalent radius
        newatom_rad = etab.GetVdwRad(newatom_atomic_num) 
    elif(bonded == "no"):
        #newatom_rad = etab.GetVdwRad(newatom_atomic_num) # get van der Waals radius
       # newatom_rad = etab.GetCovalentRad(newatom_atomic_num) # get covalent radius
       # newatom_rad = etab.GetCovalentRad(1) # use hydrogen covalent radius
        newatom_rad = etab.GetVdwRad(1)
    newatom_xyzr = [0,0,0,newatom_rad]
    
    # set box dimensions based on molecule size    
    box = Box(xyzr,0.5,10)
    
    # create accessibilty map
    
    probe_radius = newatom_rad 
    acc_map = box.createAtomAccessibiltyMap(probe_radius)
    print "addAtomToMol: Finished creating atom accessibility map."
    xc = xyzr[near_atomnum-1][0]
    yc = xyzr[near_atomnum-1][1]
    zc = xyzr[near_atomnum-1][2]
    rc = xyzr[near_atomnum-1][3]
    

#    # find the box edges, the atom is close to
#    dist1 = numpy.abs(xc-box.xleft)
#    dist2 = numpy.abs(xc-box.xright)
#    xmin=box.xleft
#    xo = -1.0
#    if(dist1>dist2):
#        xmin=box.xright
#        xo = 1.0
#    dist1 = numpy.abs(yc-box.yleft)
#    dist2 = numpy.abs(yc-box.yright)
#    ymin=box.yleft
#    yo = -1.0
#    if(dist1>dist2):
#        ymin=box.yright
#        yo = 1.0
#    dist1 = numpy.abs(zc-box.zleft)
#    dist2 = numpy.abs(zc-box.zright)
#    zmin=box.zleft
#    zo = -1.0
#    if(dist1>dist2):
#        zmin=box.zright
#        zo = 1.0
#    rdist = rc + probe_radius -0.1
    rdist = rc + probe_radius
    print "rdist = ", rdist
    
    # randomly select a position for the new atom near the specified atom of the
    # molecule

    found = 0
    delr = 0.01
    ntrials=5000
    ndist = 50    
    count = 0

    print "Selecting a position for the new atom near the atom located at"
    print xc, yc, zc
    
    while (found==0 and count < ndist):
        print "count #", count
        rdist = rdist+delr
        ntrials = 5000
        trial = 0
        
        x1 = xc - rdist
        x2 = xc + rdist
        if (x1<box.xleft):
            x1 = box.xleft
        
        if(x2>box.xright):
            x2 = box.xright
        
        y1 = yc - rdist
        y2 = yc + rdist
        if (y1<box.yleft):
            y1 = box.yleft
        
        if(y2>box.yright):
            y2 = box.yright
    
        z1 = zc - rdist
        z2 = zc + rdist
        if (z1<box.zleft):
            z1 = box.zleft
        
        if(z2>box.zright):
            z2 = box.zright
                
        print "x1, y1, z1: ", x1, y1, z1
        print "x2, y2, z2: ", x2, y2, z2
        while (trial < ntrials and found==0):
            print "trial #: ", numpy.str(trial)
            rand_x = x1+numpy.random.ranf()*(x2-x1)
            rand_y = y1+numpy.random.ranf()*(y2-y1)
            rand_z = z1+numpy.random.ranf()*(z2-z1)
            
            ix = numpy.int((rand_x-box.xleft)/box.h) 
            iy = numpy.int((rand_y-box.yleft)/box.h) 
            iz = numpy.int((rand_z-box.zleft)/box.h) 
            
            if(acc_map[ix,iy,iz]==0):
                 newatom_xyzr[0] = rand_x
                 newatom_xyzr[1] = rand_y
                 newatom_xyzr[2] = rand_z
                 found = 1
                 print "Found a site, ", numpy.str(rand_x), numpy.str(rand_y),
                 print numpy.str(rand_z)
                 rd = numpy.sqrt((rand_x-xc)*(rand_x-xc)+(rand_y-yc)*(rand_y-yc)+
                             (rand_z-zc)*(rand_z-zc))
                 print "inter-atom distance (A): ",rd
                
            else:
                trial = trial + 1
                found = 0
            
                
        count = count + 1
    
    iatom = mol.atoms[near_atomnum-1]
    iatom_atomicnum = iatom.atomicnum
    
    
    
    if(found==0):
        msg = "[geometry.py(addAtomToMol)]: Could not find a site to place "
        msg + add_element
        logger.info(msg)
        
        sys.exit(msg)
    else:
       
        msg = "added " + add_element + " at a distance = "
        msg = msg + numpy.str(rd) + " from atom " + etab.GetName(iatom_atomicnum)
        msg = msg + " (atom number = " + str(near_atomnum) + ")"
    
        logger.info(msg)
    
        
        newatom = mol.OBMol.NewAtom()
        newatom.SetAtomicNum(newatom_atomic_num)
        newatom.SetVector(rand_x,rand_y,rand_z)
        
        xyzr.append(newatom_xyzr)
        
        if(covalent_bond_type == "single" and bonded == "yes"):
            new_atomnum = len(mol.atoms)
            
            mol.OBMol.AddBond(near_atomnum,new_atomnum,1)
            
            # good for single bonds. For double and triple bonds, we have to
            # the appropriate covalent radius. For our purposes, we only need 
            # single bond (for protonation).
            
    return(mol,xyzr)
    
        
        
#Function to remove atom from molecule
def removeAtomFromMol(mol,atom_index):
    """
    This function removes an atom from the pybel moleule object, and saves it
    in .mol file format.
    INPUT:
        mol:  pybel molecule object 
        inchi_key: inChI key
        atom_index: index number of the atom in the atom list of the molecule
      
    
    OUTPUT:
        pybelmol: the new pybel molecule object
        total_chg: total charge on the new molecule
        mol3Dfile:  the .mol file of the adduct
    """
    
    # mol should contain the atom coordinates
    rematom = mol.OBMol.GetAtomById(atom_index) # get the pointer to the atom to be removed
    #atomic_num = rematom.GetAtomicNum() # get the atomic number of the atom to be removed

    #etab = openbabel.OBElementTable()    
   # element_symbol = etab.GetSymbol(atomic_num) # get the element symbol
    
    mol.OBMol.DeleteAtom(rematom) # delete the atom
    pybelmol = pybel.Molecule(mol)  # this step might be unnecessary, as the 
                                    # mol object is assumed to be pybel molecule object
    
   
    # calculate total charge on the molecule based on mmff94 forcefield
#    chgs = pybelmol.calccharges(model="mmff94")
#    total_chg = 0
#    for chg in chgs:
#       total_chg  = total_chg + chg
            
    total_chg = "NA"
    return(pybelmol,total_chg)
    
# Class for representing box 
class Box:
    """
    This class represents the dimensions based on the molecule geometry
    """
    def __init__(self,xyzr,h = 0.5,ext = 4):
        self.xmin = []
        self.xmax = []
        self.ymin = []
        self.ymax = []
        self.zmin = []
        self.zmax = []
        self.h = h
        self.ext = ext
        self.nx = []
        self.ny= []
        self.nz = []
        self.x = []
        self.y = []
        self.z = []
        self.xyzr = xyzr
    
    
        """
        This function sets the box dimensions based on the molecule geometry.
        INPUT:
            h: grid spacing in x,y, and z directions of the mesh on to
                         which the molecule geometry is mapped.
            xyzr: array of x,y,z, and radius values for each atom in the molecule
            ext:  distance of the box edge from the atom closest to the edge.
     
        OUTPUT:
            pybelmol: the new pybel molecule object
            total_chg: total charge on the new molecule
            mol3Dfile:  the .mol file of the adduct
        """
            # Determine the edge positions of the box along the x-axis
        [xleft,xright] = self.getBoxEdges(0)
    
        # Determine the edge positions of the box along the x-axis
        [yleft,yright] = self.getBoxEdges(1)
        
        # Determine the edge positions of the box along the x-axis
        [zleft,zright] = self.getBoxEdges(2)
        
        
        nx = numpy.int(numpy.round((xright - xleft)/self.h) + 1)
        ny = numpy.int(numpy.round((yright - yleft)/self.h) + 1)
        nz = numpy.int(numpy.round((zright - zleft)/self.h) + 1)
        
        print "nx =", nx, ",ny =", ny, ",nz =", nz
        
        # Correct the values of edge distance based on nx, ny, and nz values.
        xright = xleft + (nx-1)*self.h
        yright = yleft + (ny-1)*self.h
        zright = zleft + (nz-1)*self.h
        
        self.xleft = xleft
        self.xright = xright
        self.yleft = yleft
        self.yright = yright
        self.zleft = zleft
        self.zright = zright
        
        self.nx = nx
        self.ny = ny
        self.nz = nz
        
        for i in range(0,nx):
            xi = self.xleft + i*self.h
            self.x.append(xi)
        
        for i in range(0,ny):
            yi = self.yleft + i*self.h
            self.y.append(yi)
            
        for i in range(0,nz):
            zi = self.zleft + i*self.h
            self.z.append(zi)
    
                    
        print "xleft = " + numpy.str(xleft) + ", and xright = " + numpy.str(xright)
        print "yleft = " + numpy.str(yleft) + ", and yright = " + numpy.str(yright)
        print "zleft = " + numpy.str(zleft) + ", and zright = " + numpy.str(zright)
    
    
    # Function to determine the box length dimension along an axis direction
    def getBoxEdges(self,index):
          # Determine the edge positions of the box along the axis
        natoms = len(self.xyzr)    
        patom_min = self.xyzr[0][index] - self.xyzr[0][3] - self.ext
        patom_max = self.xyzr[0][index] + self.xyzr[0][3] + self.ext
        
        for i in range(1,natoms): # i:1, 2, ...,natoms - 1
            tmp = self.xyzr[i][index] + self.xyzr[i][3] + self.ext # tmp = p + rad + ext
            if(patom_max < tmp):
                patom_max = tmp
            tmp = self.xyzr[i][index] - self.xyzr[i][3] - self.ext # tmp = p - rad - ext
            if(patom_min > tmp):
                patom_min = tmp
        
        left_ex = patom_min/self.h
        right_ex = patom_max/self.h
        
        pleft = numpy.floor(left_ex)*self.h - self.ext
        pright = numpy.ceil(right_ex)*self.h + self.ext
        return(pleft,pright)
    
    
    # create atom accessibility map
    def createAtomAccessibiltyMap(self,prob_rad):
        """
        This function creates the molecule surface accessibility map of the new 
        atom to be added in the molecule.
        INPUT:
            prob_rad: radius of the atom to be added
     
        OUTPUT:
            pybelmol: the new pybel molecule object
            total_chg: total charge on the new molecule
            mol3Dfile:  the .mol file of the adduct
        """
        
        # initialize the values of the accessibilty map to zero
        acc_map = numpy.zeros((self.nx,self.ny,self.nz))
        
        natoms = len(self.xyzr)
        
#        
#        Create the  accessibility map by
#        rolling a spherical particle of radius (= probe radius) on the surface of each atom
#        to determine the domain accessible/in-accessible to the particle.
 
        for i in range(0,natoms):
            atom_radius = self.xyzr[i][3]
            atom_x = self.xyzr[i][0]
            atom_y = self.xyzr[i][1]
            atom_z = self.xyzr[i][2]
            
            r_atom = atom_radius + prob_rad
            
            x1 = numpy.floor((atom_x-self.xleft-r_atom)/self.h)
            x2 = numpy.ceil((atom_x-self.xleft+r_atom)/self.h)
            ix1 = numpy.int(x1)
            ix2 = numpy.int(x2)
            
            y1 = numpy.floor((atom_y-self.yleft-r_atom)/self.h)
            y2 = numpy.ceil((atom_y-self.yleft+r_atom)/self.h)
            iy1 = numpy.int(y1)
            iy2 = numpy.int(y2)
            
            z1 = numpy.floor((atom_z-self.zleft-r_atom)/self.h)
            z2 = numpy.ceil((atom_z-self.zleft+r_atom)/self.h)
            iz1 = numpy.int(z1)
            iz2 = numpy.int(z2)
            
            for ix in range(ix1,ix2+1):
                for iy in range(iy1,iy2+1):
                    for iz in range(iz1,iz2+1):
                        
                        xdist = self.x[ix]-atom_x
                        ydist = self.y[iy]-atom_y
                        zdist = self.z[iz]-atom_z
                        rdist = numpy.sqrt(xdist*xdist+ydist*ydist+zdist*zdist)
                        
                        if(rdist < r_atom):
                            #print numpy.str(ix), numpy.str(iy), numpy.str(iz)
                            acc_map[ix,iy,iz] = 1.0
                        
            
        return(acc_map)


# Function to convert from inchi to mol file format using openbabel library 
# instead of pybel
def inchi2molOB(inchi_str,inchi_key,file_prefix,filedir,forcefield):
    obConversion = openbabel.OBConversion()

    obConversion.SetInAndOutFormats("inchi","mol")
    
    mol = openbabel.OBMol()
    obConversion.ReadString(mol,inchi_str)
    mol.AddHydrogens()
    
    strs = inchi_key.split("InChIKey=")
    molfile =filedir + file_prefix + strs[1] + ".mol"
    obConversion.WriteFile(mol,molfile)
   
