# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 11:21:44 2016

@author: yesi172
"""

import re
import os
import glob
import sys

import inchi2input
import exporters
import loggers

def IsFloat(sValue):
    '''Function to check a given value has in float type or not'''
    try:
        float(sValue)
        return True
    except ValueError:
        return False
    
def ReplaceText(sFileName, str1, str2):
    '''Function to replace a string to another one in a given file'''
    with open(sFileName, 'r') as file:
       lines = file.readlines()
    with open(sFileName, 'w') as file:
       for line in lines:
           file.write(re.sub(str1, str2, line))
           
def InsertText(sFileName, str1, str2):
    '''Function to insert a string next to a specified one in a given file'''
    with open(sFileName, 'r') as file:
        lines = file.readlines()
    with open(sFileName, 'w') as file:
        for line in lines:
            if str1 in line: 
                # line has new line in end, remove it after adding str2 
                # add newline
                file.write(line[:-1] + '\t' + str2 + '\n') 
            else:
                file.write(line)

def DiscardCharacter(str1, str2):
    '''Function to remove a character from a list'''
    return str1.replace(str2, "")
                      
def CheckText(sFileName, sBlock):
    
    '''Function to check if a table exists in a given file'''
    
    bCIsValid = False #does not exist

    if len(sBlock.split('\n')) >= 2: 
        sBlock1 = sBlock.split('\n')[0]
        sBlock2 = sBlock.split('\n')[1]   
        with open(sFileName, 'r') as f:    
            lines = f.readlines() 
            for i, line in enumerate(lines):
                if sBlock1 in line:
                    if sBlock2 in lines[i+1]:
                        bCIsValid = True  
    else:
        with open(sFileName, 'r') as f:    
            lines = f.readlines() 
            for i, line in enumerate(lines):
                if sBlock in line:
                    bCIsValid = True #Exist   
                    
    return  bCIsValid  

def hasNumbers(inputString):
    
    return bool(re.search(r'\d', inputString))
    
def detect_walltime_format(sWallTime):
    
    bCisValid = True
    
    if not ":" in sWallTime:
        bCisValid = False
    else:
        lsWallTime = sWallTime.split(":")
        for i in range(0, len(lsWallTime)):
            for j in range(0, len(str(lsWallTime[i]))):
                if hasNumbers(lsWallTime[i][j]) == False:
                    bCisValid = False
                
    return bCisValid
        
def detectformat(cell):
    
    if "InChI" == cell[0:5]:
        return "InChI"
    elif ".xyz" == cell[-4:]:
        return ".xyz"
    else: 
        return None
        
def method(tParameters):
    '''
    Function to keep the methods used at a given run together, which 
    forms the file names
    '''
    return str(tParameters[0].upper() + '/' + tParameters[1].upper() + '//' 
            + tParameters[2].upper() + '/' + tParameters[3].upper() + '/' 
            + tParameters[4].upper() + '/' + tParameters[5].upper())
    
def dicTables(sBaseName, sBlockName): 
    
    '''
    Function to keep data in a dictionary for a given table. 
    sFileName is the file which consists of tables.  Each table is 
    represented by a sBlockName (i.e. Experimental Shifts, 
    Isotropic Shieldings, ..).  Each blocks consists of a series of 
    nucleus and their corresponding data.
    '''    
    
    dictTables=dict()
    dictTables[sBlockName] = list()
    dictTableNames = [sBlockName]
    table = dictTableNames[0]
    
    with open(sBaseName + "_3D.mol", 'r') as f1:  
        lines = f1.readlines()
    
    try:
        for table in dictTables:        
            lineCounter=0
            while lineCounter < len(lines):
                line = lines[lineCounter]
                if ((len(table.split("\n")) >=2 and 
                    table.split("\n")[1] in lines[lineCounter+1] and 
                    table.split("\n")[0] in lines[lineCounter]) or 
                    (len(table.split("\n")) ==1 and 
                    table.split("\n")[0] in lines[lineCounter])):
                    while lineCounter < len(lines):
                        line=lines[lineCounter]
                        ls = line.split()
                        while len(ls) != 0 and ls[0].isdigit() == True:                 
                            dictTables[table].append(line[:-1].split())
                            lineCounter = lineCounter+1
                            if lineCounter == len(lines):
                                return dictTables[sBlockName]
                            if len(lines[lineCounter].split()) == 0:
                                return dictTables[sBlockName]
                            line=lines[lineCounter]
                            ls = line.split()
                        lineCounter = lineCounter+1
                lineCounter = lineCounter+1
    except:
        print ("\nWARNING! Table could not be found in the file: " 
                + sBaseName + "_3D.mol for\n" + sBlockName)
        
    return dictTables[sBlockName]
     
def CheckExistLocal(sFileName):
    
    '''
    This funtion checks if a given file exists in local computer.
    sFileName includes the file directory.
    '''
   
    bCIsValid = True #exist

    if not "*" in sFileName:        
        if not os.path.exists(sFileName):
            bCIsValid = False # does not exist
    else:    
        FileList = glob.glob(sFileName)
        if len(FileList) == 0: 
            bCIsValid = False # does not exist
    
    return bCIsValid   

def CheckSeriesCompleted(sLocalDir, tParameters, extension, sheet):
    
    '''
    This function takes a list of molecules in an Excel sheet and checks 
    if the output file of each molecule exists in local computer.
    '''
    
    count = 0    
    bCIsValid = False    
    
    for irow in xrange(0,sheet.nrows):
        
        cell_input = sheet.cell_value(irow, 0).encode('utf-8')
        inputtype = detectformat(cell_input)
        if inputtype == "InChI":
            InChI, InChIKey = inchi2input.InchIs(sheet, irow, 0)   
            sBaseName = InChIKey 
        elif inputtype == ".xyz":
            sBaseName = cell_input.split(".xyz")[0] 
        
        if CheckExistLocal(sLocalDir + "/" + sBaseName + "_" + 
                            str(method(tParameters).replace("/","_")) + 
                            ".output*") == True:
            count += 1
            
    if count == sheet.nrows:
        bCIsValid = True

    return bCIsValid

def CheckIfAllOutputsCompleted(sheet1, sheet2, sLocalDir):
    
    '''
    This function takes a list of molecules in an Excel sheet as sheet1. 
    It takes a list of methods in an Excel sheet as sheet3. 
    It checks if the output files of molecules for each method exist 
    in local computer or not.
    '''
    
    counter = 0      
    for i in range(1, sheet2.nrows):
        
        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i) 
        
        if CheckSeriesCompleted(sLocalDir, tParameters, ".output*", sheet1) == True:
            counter += 1 
         
    return counter
            
def CheckOutputCompletedLocal(sFileName):   
    
    '''
    This function checks if the simulation is completed for a molecule 
    by checking its output file.
    '''    

    if CheckText(sFileName, "Total times  cpu:") == True: 
        return True
    else:
        return False

def RemoveInCompleteOutputFileLocal(sFileName):   

    '''
    This function removes the given file in local computer if it is incomplete 
    '''    
    if CheckExistLocal(sFileName) == True:

        if CheckOutputCompletedLocal(sFileName) == False: 

            os.remove(sFileName)

def RemoveInCompleteOutputFilesLocal(sheet1, sheet2, sLocalDir):   
    
    '''
    This function removes incomplete output files in local computer
    '''    
    
    for i in range(1, sheet2.nrows):
        
        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i)
        
        inputtype = detectformat(tParameters[6])
        if inputtype == "InChI":
            refInChI = tParameters[6]        
            sRefBaseName = inchi2input.inchi2inchikey(refInChI)
        elif inputtype == ".xyz":
            sRefBaseName = tParameters[6].split(".xyz")[0]
        
        sFileName = (sLocalDir + "\\" + sRefBaseName + "_" + 
                                str(method(tParameters).replace("/","_")) + 
                                ".output")
        RemoveInCompleteOutputFileLocal(sFileName)

        for irow in xrange(0, sheet1.nrows):
            
            cell_input = sheet1.cell_value(irow, 0).encode('utf-8')
            inputtype = detectformat(cell_input)
            if inputtype == "InChI":
                InChI, InChIKey = inchi2input.InchIs(sheet1, irow, 0)   
                sBaseName = InChIKey 
            elif inputtype == ".xyz":
                sBaseName = cell_input.split(".xyz")[0] 
            
            sFileName = (sLocalDir + "\\" + sBaseName + "_" + 
                                str(method(tParameters).replace("/","_")) + 
                                ".output")            
            RemoveInCompleteOutputFileLocal(sFileName)
       
def CheckMustParser(sLocalDir, sFileName1, sFileName2, sHost, sServerDir):
    
    '''
    This function raises error messages for the must parsers if not 
    given by the user.
    '''
    
    #if sLocalDir is None: 
    #    sys.exit("\nFatal Error: File directory in local computer was not " 
    #            + "specified!\n\nPlease enter full path of working directory "
    #            + "in Local Computer. "
    #            + "(Ex: -dir C:/Users/UserName/Desktop/test/)\n\n"
    #            + "For more information please type python isicle.py --help")
    if sFileName1 is None: 
        sys.exit("\nFatal Error: Excel file with InchI set was not specified!"
                + "\n\nPlease enter name of Excel file in Local Computer. "
                + "(Ex: -molcs Molecule_Set.xlsx)\n\n"
                + "For more information please type python isicle.py --help")
    if sFileName2 is None: 
        sys.exit("\nFatal Error: Excel file for DFT Methods was not specified!"
                + "\n\nPlease enter name of Excel file in Local Computer. "
                + "(Ex: - metd Methods.xlsx)\n\n"
                + "For more information please type python isicle.py --help")
    if sHost is None: 
        sys.exit("\nFatal Error: Host name in Cascade was not specified!\n\n"
                + "Please enter Account name in Cascade. "
                + "(Ex: -host NetworkID@cascade.emsl.pnl.gov)\n\n"
                + "For more information please type python isicle.py --help")
    if sServerDir is None: 
        sys.exit("\nFatal Error: Destination path was not specified!\n\n"
                + "Please enter full path of working directory in Cascade. "
                + "(Ex: -dest /dtemp/NetworkID/)\n\n"
                + "For more information please type python isicle.py --help")
        
def CheckReferenceMolecule(sLocalDir, sheet2):
    
    '''
    This function checks a possible future error regarding the 
    reference molecule 
    '''
    
    for i in range(1, sheet2.nrows):
            
        tOriginalParameters, tParameters  = exporters.ExtractInputs(sheet2, i)
    
        inputtype = detectformat(tParameters[6])
        
        if inputtype == None:
            sys.exit("ERROR in Method no: " + str(i) + method(tOriginalParameters)
                    + "format of reference compound cannot be deteced!")
        
        if inputtype == "InChI":
            refInChI = tParameters[6]        
            refInChIKey = inchi2input.inchi2inchikey(refInChI)
            sysErrorMsgInChI(sLocalDir, refInChI, refInChIKey, tParameters[4])
        elif inputtype == ".xyz":
            sRef_filename = tParameters[6]
            sysErrorMsgXYZ(sLocalDir, sRef_filename, tParameters[4])
        
def sysErrorMsgInChI(sLocalDir, InChI, InChIKey, nuclei):
    
    '''
    This function creates the 3D mol file for reference molecule 
    if does not exist.  It also gives an error message if the nucleus 
    or nuclei have been looking for are not found in reference compound. 
    '''
    
    if CheckExistLocal(sLocalDir + "/" + InChIKey 
                      + "_3D.mol") == False: # if not exist
       # Create 3Dmol file
       inchi2input.inchi2triDmol(InChI, InChIKey, sLocalDir, 'mmff94')
       inchi2input.inchi2XYZFile(InChI, InChIKey, sLocalDir, 'mmff94')
       loggers.AddExpState(InChIKey)
    elif CheckText(sLocalDir + "/" + InChIKey + "_3D.mol", 
                   "Experimental Shifts") == False:
            loggers.AddExpState(InChIKey)
            
    if len(exporters.NucleusNoListXYZ(InChIKey, nuclei)) == 0: 
        sys.exit("Error msg: Reference compound does not have "
                + "the specified nucleus")
    
    for i in range(0, len(list(nuclei))):  
        if not list(nuclei)[i] in exporters.ReadXYZFile(InChIKey + ".xyz")[1]:
            sys.exit("Error msg: Reference compound does not have "
                    + "the specified nucleus")

def sysErrorMsgXYZ(sLocalDir, sRef_filename, nuclei):
    
    '''
    This function creates the 3D mol file for reference molecule 
    if does not exist.  It also gives an error message if the nucleus 
    or nuclei that have been looking for are not found in reference compound. 
    '''
    
    sFileName = sLocalDir + "/" + sRef_filename.split(".xyz")[0] + "_3D.mol"    
    if CheckExistLocal(sFileName) == False: # if not exist 
        loggers.AddExpState(sRef_filename.split(".xyz")[0])
     
    sFileName = sLocalDir + "/" + sRef_filename
    if len(exporters.ReadXYZFile(sFileName)[1]) == 0: 
        sys.exit("Error msg: Reference compound does not have "
                + "the specified nucleus")   

    for i in range(0, len(list(nuclei))):
        if not list(nuclei)[i] in exporters.ReadXYZFile(sFileName)[1]:
            sys.exit("Error msg: Reference compound does not have "
                    + "the specified nucleus")   
