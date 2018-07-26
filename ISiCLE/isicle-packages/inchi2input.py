# -*- coding: utf-8 -*-
# Author: Yasemin Yesiltepe

import os
import shutil

import pybel
import openbabel

import tools
import exporters
import loggers

def inchi2molOB(InChI, InChIKey, sLocalDir):

    '''
    Function to convert from inchi to mol file format using Openbabel 
    library instead of pybel
    '''

    mol = openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("inchi","mol")
    obConversion.ReadString(mol,InChI)    
    mol.AddHydrogens()
    # File name for mol file
    sFileName = sLocalDir + "\\" + InChIKey + ".mol" 

    # outMDL = obConversion.WriteString(mol)
    # obConversion.SetOptions("gen3D", obConversion.Option_type.OUTOPTIONS)
    obConversion.WriteFile(mol, sFileName)

def smi2molOB(Smiles, InChIKey, sLocalDir):

    '''
    Function to convert from smiles to mol file format using Openbabel 
    library instead of pybel
    '''

    mol = openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi","mdl")
    obConversion.ReadString(mol, Smiles)
    mol.AddHydrogens()
    # File name for mol file
    sFileName = sLocalDir + "\\" + InChIKey + ".mol" 
    # outMDL = obConversion.WriteString(mol)

    obConversion.WriteFile(mol, sFileName)

def inchi2triDmolOB(InChI, InChIKey, sLocalDir):

    '''
    Function to convert from inchi to 3D mol file format using openbabel 
    library instead of pybel
    '''

    mol = openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("inchi","mol")
    obConversion.ReadString(mol, InChI)
    mol.AddHydrogens()
    # mol.Center()
    
    # File name for mol file
    sFileNameInput = sLocalDir + "\\" + InChIKey + ".mol" 
    # File name for 3D mol file
    sFileNameOutput = sLocalDir + "\\" + InChIKey + "_3D.mol" 
    
    #outMDL = obConversion.WriteString(mol)
    obConversion.WriteFile(mol, sFileNameInput)    

    os.system("obabel " + sFileNameInput + " -O " 
             + sFileNameOutput + " --gen3d -c")

def inchi2xyzOB(InChI, InChIKey, sLocalDir):

    '''
    Function to convert from inchi to xyz file format using openbabel 
    library instead of pybel
    '''

    mol = openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("inchi","mol")
    obConversion.ReadString(mol, InChI)
    mol.AddHydrogens()
    #mol.Center()
    
    # File name for mol file
    sFileNameInput = sLocalDir + "\\" + InChIKey + ".mol" 
    # File name for xyz file
    sFileNameOutput = sLocalDir + "\\" + InChIKey + ".xyz" 
    
    #outMDL = obConversion.WriteString(mol)
    obConversion.WriteFile(mol, sFileNameInput)

    os.system("obabel " + sFileNameInput + " -O " 
             + sFileNameOutput + " --gen3d -c")

def smi2xyzOB(Smiles, InChIKey, sLocalDir):

    '''
    Function to convert from smiles to xyz file format using 
    openbabel library instead of pybel
    '''

    mol = openbabel.OBMol()
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi","mol")
    obConversion.ReadString(mol, Smiles)
    mol.AddHydrogens()
    # mol.Center()

    # File name for mol file
    sFileNameInput = sLocalDir + "\\" + InChIKey + ".mol" 
    # File name for xyz file
    sFileNameOutput = sLocalDir + "\\" + InChIKey + ".xyz" 
    
    # outMDL = obConversion.WriteString(mol)
    obConversion.WriteFile(mol, sFileNameInput)

    os.system("obabel " + sFileNameInput + " -O " 
             + sFileNameOutput + " --gen3d -c")

def minimizeOB(sFileName, ffield = "mmff", iSteps = 300, fTolerance = 1e-5):

    '''
    optimize the geometry, minimize the energy for a molecule
    Synopsis: obminimize [OPTIONS] filename
    Options: If no filename is given, obminimize will give all options 
    including the available forcefields.

    -n steps
    Specify the maximum number of steps (default=2500)
    -cg
    Use conjugate gradients algorithm (default)
    -sd
    Use steepest descent algorithm
    -c criteria
    Set convergence criteria (default=1e-6)
    -ff forcefield
    Select the forcefield
    '''

    # Minimize the energy for the molecule(s) in file test.mol2
    os.system("obminimize " + sFileName)

    # Minimize the energy for the molecule(s) in file test.mol2 using 
    # the given forcefield:
    # os.system("obminimize -ff " + ffield + " " + sFileName)

    # Minimize the energy for the molecule(s) in file test.mol2 by taking 
    # at most 300 geometry optimization steps
    # os.system("obminimize -n " + iSteps + " " + sFileName)

    # Minimize the energy for the molecule(s) in file test.mol2 using 
    # the steepest descent algorithm and convergence criteria 1e-5:
    # os.system("obminimize -sd -c " + fTolerance + " " + sFileName)

def check_substructure(InChI, sSmarts):
    
    molecule = pybel.readstring("inchi", InChI)
    smarts = pybel.Smarts(sSmarts) 
    
    return smarts.findall(molecule)

def inchiList2sdf(InchIList):
    
    molecules = [pybel.readstring("inchi", x) for x in InchIList]
    output = pybel.Outputfile("sdf", "outputfile.sdf", overwrite=False)
    for i in xrange(0, len(molecules)):
        output.write(molecules[i])
    output.close()
        
def xyz2mol(sBaseName):
    
    sFileNameInput = sBaseName + ".xyz"
    sFileNameOutput = sBaseName + ".mol"
    
    os.system("babel -ixyz " + sFileNameInput + " -omol " + sFileNameOutput)

def smi2inchi(Smiles):
    
    '''Function to convert Smiles to InChI code'''
    
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("smi", "inchi")

    mol = openbabel.OBMol()
    conv.ReadString(mol, Smiles)
    conv.SetOptions("I", conv.OUTOPTIONS)
    InChI = conv.WriteString(mol).rstrip()

    return InChI.encode('utf-8')

def bIsInRing(sFileNameInput, sBaseName, sAtomNumber):
    
    iAtomNumber = int(sAtomNumber) - 1
    
    for molecule in pybel.readfile("sdf", sFileNameInput):
        data = molecule.data #Type: pybel.MoleculeData
        if data.has_key("INCHI_KEY") and sBaseName == data['INCHI_KEY']:
            return molecule.atoms[iAtomNumber].OBAtom.IsInRing()

def AtomType(sFileNameInput, sBaseName, sAtomNumber):

    iAtomNumber = int(sAtomNumber) - 1
    
    for molecule in pybel.readfile("sdf", sFileNameInput):
        data = molecule.data #Type: pybel.MoleculeData
        if data.has_key("INCHI_KEY") and sBaseName == data['INCHI_KEY']:
            return molecule.atoms[iAtomNumber].type

def AtomData(sFileNameInput, sBaseName, sAtomNumber):
    
    iAtomNumber = int(sAtomNumber) - 1
    
    for molecule in pybel.readfile("sdf", sFileNameInput):
        data = molecule.data #Type: pybel.MoleculeData
        if data.has_key("INCHI_KEY") and sBaseName == data['INCHI_KEY']:
            return sAtomNumber + " ", molecule.atoms[iAtomNumber].idx, molecule.atoms[iAtomNumber].type, molecule.atoms[iAtomNumber].OBAtom.IsInRing()

def ReadXYZFile(sFileName):
    '''
    It reads the given xyz file.  Filename includes extension of .xyz. 
    Ex: UHOVQNZJYSORNB-UHFFFAOYSA-N.xyz  It skips the first two lines. 
    First line is the total number of atoms in the file. Second line 
    says the program which the file was generated by.  
    
    Format of Input file:  
    3
    generated by VMD
    O         2.043000       -1.222000       -0.125000
    H         1.599000       -1.723000        0.610000
    H         2.933000       -1.531000       -0.038000
    
    Output file is a list of atoms and xyz coordinates in a string type. 
    Format of Output: 
    [['O', '2.043000', '-1.222000', '-0.125000'],
     ['H', '1.599000', '-1.723000', '0.610000'],
     ['H', '2.933000', '-1.531000', '-0.038000']]
    '''
    lXYZFile = list()
     
    with open(sFileName, 'r') as f:    
            lines = f.readlines() 
            for i, line in enumerate(lines):
                if i >= 2: #it skips the first two lines (first line is the number of atoms in the file    generated by VMD) 
                    lXYZFile.append(line.split()) 
    
    if lXYZFile[-1] == []: del lXYZFile[-1]
        
    lsXYZCoorMolecule = [(lXYZFile[row][1:]) for row in xrange(0,len(lXYZFile))]
    lsXYZAtomMolecule = [(lXYZFile[row][0]) for row in xrange(0,len(lXYZFile))]
    lfXYZCoorMolecule = [float(lsXYZCoorMolecule[row][column]) for column in xrange(0,len(lsXYZCoorMolecule[row])) for row in xrange(0,len(lXYZFile))]
    lfXYZCoorMolecule = numpy.reshape(lfXYZCoorMolecule, (len(lsXYZCoorMolecule), -1), order='F')
                    
    return lXYZFile, lsXYZAtomMolecule, lsXYZCoorMolecule, lfXYZCoorMolecule #(type: list)  

def ReadMolFile(sFileName):
    '''
    It reads the given mol file. Filename includes extension of .mol. 
    Ex: UHOVQNZJYSORNB-UHFFFAOYSA-N.mol
    It skips the first two lines. 
    First line is the total number of atoms in the file. 
    Second line says the program which the file was generated by.  
    
    Format of Input file:      
    \n
     OpenBabel09211601333D
    \n
      3  2  0  0  0  0  0  0  0  0999 V2000
        1.0500   -0.0003    0.0866 C   0  0  0  0  0  0  0  0  0  0  0  0
       -0.1470   -0.0003    0.0866 O   0  0  0  0  0  0  0  0  0  0  0  0
        2.2470   -0.0003    0.0866 O   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  2  0  0  0  0
      1  3  2  0  0  0  0
    M  END
    
    Output file is a list of atoms and xyz coordinates in a string type. 
    Format of Output: 
    [['C', '1.0500', '-0.0003', '0.0866'],
     ['O', '-0.1470', '-0.0003', '0.0866'],
     ['O', '2.2470', '-0.0003', '0.0866']]
    '''
    #it skips the first four lines (first line is the number of atoms in the file generated by VMD) 
    lMolFile = list()
     
    with open(sFileName, 'r') as f:    
            lines = f.readlines() 
            for i, line in enumerate(lines):
                if "M  END\n" in line:
                    break
                if i >= 4 and len(line) >= 30: 
                    lMolFile.append([line.split()[3], line.split()[0], line.split()[1], line.split()[2]]) 

    lsMolCoorMolecule = [(lMolFile[row][1:]) for row in xrange(0,len(lMolFile))]
    lsMolAtomMolecule = [(lMolFile[row][0]) for row in xrange(0,len(lMolFile))]
    lfMolCoorMolecule = [float(lsMolCoorMolecule[row][column]) for column in xrange(0,len(lsMolCoorMolecule[row])) for row in xrange(0,len(lMolFile))]
    lfMolCoorMolecule = numpy.reshape(lfMolCoorMolecule, (len(lsMolCoorMolecule), -1), order='F')
                    
    return lMolFile, lsMolAtomMolecule, lsMolCoorMolecule, lfMolCoorMolecule #(type: list) 
    
def MakeCenter(lfVectors, lCenter):
    
    xyzdistance = list()
    Center = GeometricCenter(lfVectors)
    xyzdistance.append(Center[0]-lCenter[0])
    xyzdistance.append(Center[1]-lCenter[1])
    xyzdistance.append(Center[2]-lCenter[2])
    
    lfVectorsTransformed = MoveBy(lfVectors, [xyzdistance[0],xyzdistance[1],xyzdistance[2]])
        
    return lfVectorsTransformed

def vector(lpoint1, lpoint2):
    '''
    It creates a vector from two points. 
    Inputs: Two lists with 3 elements each. It takes float numbers.
    Ex: lpoint1: [2.043000, -1.222000, -0.125000], lpoint2: [1.599000, -1.723000, 0.610000]
    '''
    lvector = [lpoint2[0]-lpoint1[0],lpoint2[1]-lpoint1[1],lpoint2[2]-lpoint1[2]]
    
    return lvector #(type: list, float) 

def Rotate(lfVectors, vector):
    
    return MoveBy(lfVectors, vector)
    
def GeometricCenter(lfvectors):
    '''
    It calculates the geometric center of a given list of xyz cartesian coordinates.
    Format of lvectors: 
    [[2.043000, -1.222000, -0.125000],
     [1.599000, -1.723000, 0.610000],
     [2.933000, -1.531000, -0.038000]]
    Output : [2.191666666666667, -1.492, 0.149]
    '''
    x = 0.
    y = 0.
    z = 0.    
    for i in xrange(0, len(lfvectors)):
        x = x + lfvectors[i][0]
    x = x / len(lfvectors)
    for i in xrange(0, len(lfvectors)):
        y = y + lfvectors[i][1]
    y = y / len(lfvectors)
    for i in xrange(0, len(lfvectors)):
        z = z + lfvectors[i][2]
    z = z / len(lfvectors)
    
    return [x, y, z] #(type: list, float) 

def MoveBy(lfVectors, xyzdistance):
    
    lfVectorsTransformed = list()
    for row in xrange(0, len(lfVectors)):
        x = lfVectors[row][0] - xyzdistance[0]
        y = lfVectors[row][1] - xyzdistance[1]
        z = lfVectors[row][2] - xyzdistance[2]
        lfVectorsTransformed.append([x, y, z])
    
    return lfVectorsTransformed

def Align(lfVectors1, lfVectors2):
    
    #align first atoms     
    rotationvector = vector(lfVectors1[0], lfVectors2[0])
    
    #keep first molecule constant, #rotate second molecule 
    lfVectors2 = Rotate(lfVectors2, rotationvector)
    
    return lfVectors1, lfVectors2

def rmsd(lfVectors1, lfVectors2):    
    deviation = sum(squared_distance(atomA, atomB) for 
                    (atomA, atomB) in zip(lfVectors1, lfVectors2))
    return math.sqrt(deviation / float(len(lfVectors1)))
    
def squared_distance(coordsA, coordsB):
    """Find the squared distance between two 3-tuples"""
    sqrdist = sum( (a-b)**2 for a, b in zip(coordsA, coordsB))    
    return sqrdist
    
def rmsdMolXYZ(sFileNameMol, sFileNameXYZ):

    lMolOutputs = ReadMolFile(sFileNameMol)
    lsMolAtomMolecule, lfMolCoorMolecule = lMolOutputs[1], lMolOutputs[3]
    
    lXYZOutputs = ReadXYZFile(sFileNameXYZ)
    lsXYZAtomMolecule, lfXYZCoorMolecule = lXYZOutputs[1], lXYZOutputs[3]
    
    #Shift coordinates by center of mass
    lfMolCoorMolecule = MoveBy(lfMolCoorMolecule, GeometricCenter(lfMolCoorMolecule))
    lfXYZCoorMolecule = MoveBy(lfXYZCoorMolecule, GeometricCenter(lfXYZCoorMolecule))

    #Move molecules to the origin
    lfMolCoorMolecule = MoveBy(lfMolCoorMolecule, [0, 0, 0])
    lfXYZCoorMolecule = MoveBy(lfXYZCoorMolecule, [0, 0, 0])

    #Align second molecule to first molecule
    lfMolCoorMolecule, lfXYZCoorMolecule = Align(lfMolCoorMolecule, lfXYZCoorMolecule)
    
    return rmsd(lfMolCoorMolecule, lfXYZCoorMolecule)

def AtomNumbers(InChI):
    
    molecule = pybel.readstring("inchi", InChI)
    molecule.addh()
    return len(molecule.atoms)

def MolecularWeight(InChI):

    molecule = pybel.readstring("inchi", InChI)
    molecule.addh()
    return molecule.molwt
    
def inchi2smi(InChI):

    '''Function to convert InChI code to Smiles'''
    
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("inchi", "smi")

    mol = openbabel.OBMol()
    conv.ReadString(mol, InChI)
    conv.SetOptions("I", conv.OUTOPTIONS)
    Smiles = conv.WriteString(mol).rstrip()

    return Smiles.encode('utf-8')

def inchi2mol(InChI, InChIKey, sLocalDir):

    '''Function to convert InChI code to mol file'''
    
    # File name for 2D mol file
    sFileName = sLocalDir + "\\" + InChIKey + ".mol" 

    # Read in a molecule from a string.
    # readstring(format, string)   #Required parameters: format, string
    molecule = pybel.readstring("inchi", InChI)
    # Write the molecule to a file or return a string.
    # Optional parameters: format: default is "SMI",
    # filename: default is None, overwite: default is False
    molecule.write('mol', sFileName, True)

def smi2mol(Smiles, InChIKey, sLocalDir):
    
    '''Function to convert Smiles to 2D mol file'''
    
    # File name for 2D mol file
    sFileName = sLocalDir + "\\" + InChIKey + ".mol" 

    # Read in a molecule from a string.
    # readstring(format, string)   #Required parameters: format, string
    molecule = pybel.readstring("smi", Smiles)
    # Write the molecule to a file or return a string.
    # Optional parameters: format: default is "SMI", 
    # filename: default is None, overwite: default is False
    molecule.write('mol', sFileName, True)

def smi2triDmol(Smiles, InChIKey, sLocalDir, ffield):

    '''Function to convert Smiles to 3D mol file'''
    
    # File name for 3D mol file
    sFileName = sLocalDir + "\\" + InChIKey + "_3D.mol"

    # Read in a molecule from a string.
    molecule = pybel.readstring("smi", Smiles)
    # Generate 3D coordinates.
    # Optional parameters: forcefield: default is "MMFF94",
    # steps: default is 50
    # Hydrogens are added
    molecule.make3D(forcefield=ffield, steps=50)
    # Locally optimize the coordinates.
    # Optional parameters: forcefield: default is "MMFF94",
    # steps: default is 500
    # Molecule needs to have explicit hydrogens. If not, call addh()
    molecule.localopt(forcefield=ffield, steps=500)
    # Write the molecule to a file or return a string.
    # Optional parameters: format: default is "SMI", 
    # filename: default is None, overwite: default is False
    molecule.write('mol', sFileName, True)

def inchi2triDmol(InChI, InChIKey, sLocalDir, ffield):

    '''Function to convert InChI code to 3D mol file'''
    
    # File name for 3D mol file
    sFileName = sLocalDir + "\\" + InChIKey + "_3D.mol" 

    # Read in a molecule from a string.
    molecule = pybel.readstring("inchi", InChI)
    # Generate 3D coordinates.
    # Optional parameters: forcefield: 
    # default is "MMFF94",steps: default is 50
    # Hydrogens are added
    molecule.make3D(forcefield=ffield, steps=50)
    # Locally optimize the coordinates.
    # Optional parameters: forcefield: 
    # default is "MMFF94",steps: default is 500
    # Molecule needs to have explicit hydrogens. If not, call addh()
    molecule.localopt(forcefield=ffield, steps=500)
    # Write the molecule to a file or return a string.
    # Optional parameters: format: default is "SMI", 
    # filename: default is None, overwite: default is False
    molecule.write('mol', sFileName, True)

def inchi2XYZFile(InChI, InChIKey, sLocalDir, ffield):

    """
    This function creates an .xyz file of a given InChI code.
    INPUT:
        InchI: inchi string
        InchIKey: : inchikey, output file name
        filedir: The directory where the .xyz file is saved
        ffield: force field, default:mmff94
    OUTPUT:
        mol:  pybel molecule object
        xyzfile: The .xyz file
    """

    mol = pybel.readstring("inchi", InChI)
    mol.make3D(forcefield=ffield, steps=50)
    mol.localopt(forcefield=ffield, steps=500)
    sFileName = sLocalDir + '\\' + InChIKey + ".xyz"
    mol.write("xyz", sFileName, True)

def smi2XYZFile(Smiles, InChIKey, sLocalDir, ffield):

    """
    This function creates an .xyz file of a given Smiles
    INPUT:
        Smiles: smiles string
        InchIKey: : inchikey, output file name
        filedir: The directory where the .xyz file is saved
        ffield: force field, default:mmff94
    OUTPUT:
        mol:  pybel molecule object
        xyzfile: The .xyz file
    """

    mol = pybel.readstring("smi", Smiles)
    mol.make3D(forcefield=ffield, steps=50)
    mol.localopt(forcefield=ffield, steps=500)
    sFileName = sLocalDir + '\\' + InChIKey + ".xyz"
    mol.write("xyz", sFileName, True)

def inchi2inchikey(InChI):

    '''Function to convert InChI code to InChIKey'''
    
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("inchi", "inchi")

    mol = openbabel.OBMol()
    conv.ReadString(mol, InChI)
    conv.SetOptions("K", conv.OUTOPTIONS)
    InchIKey = conv.WriteString(mol).rstrip()

    return InchIKey.encode('utf-8')
    
def MolAndXYZInputs(cell_input, sBaseName, sLocalDir):
    
    '''It generates mol and xyz files for a given InChI if not exist'''
    
    inputtype = tools.detectformat(cell_input)
    
    if inputtype == "InChI":
        InChI = cell_input
        InChIKey = inchi2inchikey(InChI) 

        # create xyz files
        inchi2XYZFile(InChI, InChIKey, sLocalDir, 'mmff94') 
        # 3D Mol file
        if (tools.CheckExistLocal(sLocalDir + "/" + InChIKey 
                                 + "_3D.mol") == False): # if not exist
            # Create 3Dmol file
            inchi2triDmol(InChI, InChIKey, sLocalDir, 'mmff94') 
            loggers.AddExpState(InChIKey)
        elif tools.CheckText(sLocalDir + "/" + InChIKey + "_3D.mol", 
                             "Experimental Shifts") == False:
            loggers.AddExpState(InChIKey)
        shutil.copyfile(sBaseName + "_3D.mol", "temp_3D.mol")
        # Create 3Dmol file
        inchi2triDmol(InChI, InChIKey, sLocalDir, 'mmff94')
        sBaseName = InChIKey 
    
    elif inputtype == ".xyz":
        sBaseName = cell_input.split(".xyz")[0]       
        if (tools.CheckExistLocal(sLocalDir + "/" + sBaseName 
                                 + "_3D.mol") == False): # if not exist
            loggers.AddExpState(sBaseName)


    #Force fields: ['UFF', 'MMFF94', MMFF94s, 'Ghemical', Gaff]

    if tools.CheckText("temp_3D.mol", 'M  END') == True:
        with open("temp_3D.mol") as f1:
            lines = f1.readlines()
            for i, line in enumerate(lines):
                if 'Experimental Shifts' in line:
                    line1 = lines[i]
                    line2 = lines[i+1]
        dictExperimentalShift = tools.dicTables("temp", "Experimental Shifts")
        with open("temp_3D.mol", 'w') as f1:
            f1.write("\n" + line1 + line2)
            for listLine in dictExperimentalShift:
                for sWord in listLine:
                    f1.write(sWord + "\t")
                f1.write("\n")

    with open("temp_3D.mol") as f:
        lines = f.readlines()
    with open(sBaseName + "_3D.mol", "a") as f1:
        for line in lines:
            f1.write(line)

def Inputs(cell_input, sBaseName, sLocalDir, 
           tOriginalParameters, tParameters, sServerDir):
    
    '''It prepares the input files to NWChem for a given molecule'''
    
    MolAndXYZInputs(cell_input, sBaseName, sLocalDir)
    # Create NWChem input files
    loggers.CreateNWCHEMFile(sBaseName, sLocalDir, sServerDir + sLocalDir.split("/")[-2], 
                             tOriginalParameters, tParameters)

def InchIs(sheet, irow, icolumn):
    
    '''It reads InChI codes of each molecule listed in the Excel file'''
    
    sInChI = sheet.cell_value(irow, icolumn).encode('utf-8')

    return sInChI, inchi2inchikey(sInChI)

def get_inputs(sLocalDir, sServerDir, cell_input, sBaseName, 
               tOriginalParameters, tParameters):
    
    inputtype = tools.detectformat(cell_input)
    
    if tools.CheckExistLocal(sLocalDir + "/" + sBaseName + ".xyz") == False: 
        # if not exist    
        if inputtype == ".xyz":
            print "ERROR! " + sBaseName + ".xyz file could not be found in the folder" 
        elif inputtype == "InChI":
            InChI = cell_input
            InChIKey = inchi2inchikey(InChI)     
            # create xyz files
            inchi2XYZFile(InChI, InChIKey, sLocalDir, 'mmff94')
   
    if tools.CheckExistLocal(sLocalDir + "/" + sBaseName + "_3D.mol") == False: 
        # if not exist
        # Create 3Dmol file
        if inputtype == "InChI":
            InChI = cell_input
            InChIKey = inchi2inchikey(InChI)     
            inchi2triDmol(InChI, InChIKey, sLocalDir, 'mmff94') 
            loggers.AddExpState(InChIKey)
        elif inputtype == ".xyz":
            loggers.AddExpState(sBaseName)
    
    # Remove incomplete output files
    if tools.CheckExistLocal(sLocalDir + "/" + sBaseName + "_" 
                    + str(tools.method(tParameters).replace("/","_")) 
                    + ".output*") == True: 
        sFileName = (sLocalDir + "/" + sBaseName + "_" 
                    + str(tools.method(tParameters).replace("/","_")) 
                    + ".output")
        tools.RemoveInCompleteOutputFileLocal(sFileName)
    
    if tools.CheckExistLocal(sLocalDir + "/" + sBaseName + "_" 
                + str(tools.method(tParameters).replace("/","_")) 
                + ".output*") == False: # if not exist
    
        if tools.CheckExistLocal(sLocalDir + "/" + sBaseName + "_" 
                    + str(tools.method(tParameters).replace("/","_")) 
                    + ".nw") == False: # if not exist
        
            loggers.CreateNWCHEMFile(sBaseName, sLocalDir,
                                                     sServerDir + sLocalDir.split("/")[-2], 
                                                     tOriginalParameters, tParameters)
 
    
def CreateAllInputs(sheet1, sheet2, sLocalDir, sServerDir,
                                sNODES, sPROCESSOR, sWalltime):
    
    '''It prepares all the input files for a given molecule set'''
    
    for i in range(1, sheet2.nrows):
        
        tOriginalParameters, tParameters  = exporters.ExtractInputs(sheet2, i)
        
        inputtype = tools.detectformat(tParameters[6])
        
        if inputtype == "InChI":
            refInChI = tParameters[6]        
            sRefBaseName = inchi2inchikey(refInChI)
        elif inputtype == ".xyz":
            sRefBaseName = tParameters[6].split(".xyz")[0]
            
        get_inputs(sLocalDir, sServerDir, tParameters[6], sRefBaseName, 
                   tOriginalParameters, tParameters)
                
        for irow in range(0, sheet1.nrows): 
            
            # number of molecules = sheet.nrows - 1
            cell_input = sheet1.cell_value(irow, 0).encode('utf-8')
            inputtype = tools.detectformat(cell_input)
            if inputtype == "InChI":
                InChI, InChIKey = InchIs(sheet1, irow, 0)   
                sBaseName = InChIKey 
            elif inputtype == ".xyz":
                sBaseName = cell_input.split(".xyz")[0]       
            
            get_inputs(sLocalDir, sServerDir, cell_input, sBaseName, 
                       tOriginalParameters, tParameters)
            print 'Input files of Molecule ' + str(irow) + ' have been successfully prepared!'
            