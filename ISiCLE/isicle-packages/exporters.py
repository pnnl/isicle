# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 11:41:53 2016

@author: yesi172
"""

import numpy
import glob
import xlrd

import inchi2input
import basicoperations
import tools

def isotropic_shieldings(sheet1, sheet2):
    
    '''
    This function writes the isotropic sheildings of all the molecules 
    given in Excel sheet
    '''
    
    for i in range(1, sheet2.nrows):
    
        tOriginalParameters, tParameters  = ExtractInputs(sheet2, i)
    
        if tools.detectformat(tParameters[6]) == "InChI":
            refInChI = tParameters[6]    
            refInChIKey = inchi2input.inchi2inchikey(refInChI)
            sRefBaseName = refInChIKey
        else:
            sRefBaseName = tParameters[6].split(".xyz")[0] 

        write_shieldings(sRefBaseName, tOriginalParameters, tParameters)     
            
            
        for irow in range(0, sheet1.nrows):
            
            cell_input = sheet1.cell_value(irow, 0).encode('utf-8')
            inputtype = tools.detectformat(cell_input)
            if inputtype == "InChI":
                InChI, InChIKey = inchi2input.InchIs(sheet1, irow, 0)   
                sBaseName = InChIKey 
            elif inputtype == ".xyz":
                sBaseName = cell_input.split(".xyz")[0] 

            write_shieldings(sBaseName, tOriginalParameters, tParameters)
            
def write_shieldings(sBaseName, tOriginalParameters, tParameters):       
    
    '''
    This function extracts the isotropic shieldings for a given 
    molecule, InChIKey, from its output file.  It writes the outputs 
    to the 3D mol files as wellas the simulation details such as 
    the number of nodes, number cores per node and cpu time.   
    '''

    if tools.CheckExistLocal(sBaseName + "_3D.mol") == False:
        
        print ("\nWARNING! Isotropic sheildings are not written "
               + "to the file \t" + sBaseName + "_3D.mol" 
               + "\nFile could not be found!")
        
        return
        
    sFileName = (sBaseName + "_"
                    + str(tools.method(tParameters).replace("/","_")))

    if tools.CheckExistLocal(sFileName + ".output") == False:
        
        print ("\nWARNING! Isotropic sheildings are not exported.\t" 
                  + sFileName + ".output" + " failed!")
        
        return

    if tools.CheckText(sFileName + ".output", "Total times") == False:

        print ("\nWARNING! Isotropic sheildings are not exported.\t" 
                  + sFileName + ".output" + " failed!")
            
        return

    if (tools.CheckText(sBaseName + "_3D.mol", "Isotropic Shielding" 
                        + "\n" + str(tools.method(tOriginalParameters))) 
                        == True):
        
        return
    
    for name in glob.glob(sFileName + ".output"):
        with open(name) as f1:
            lines = f1.readlines()  
            
    with open(sBaseName + "_3D.mol", "a") as f2:
            
        f2.write("\n\n" + "Isotropic Shielding" + "\n"
                + tools.method(tOriginalParameters) + "\n"             
                + "Number of nodes: "
                + str(Specs(sFileName)[0]) + "\t"
                + "Number of cores per node: "
                + str(Specs(sFileName)[1]) + "\t"
                + "Total CPU time: " + str(CPUtime(sFileName)) + "\t"
                + "Total time: "
                + str(Specs(sFileName)[2]) + "s"
                + "\n")
        
        try:
            lNucleusNo = read_shielding_index(sFileName + ".output")
        except:
            lNucleusNo = NucleusNoListXYZ(sBaseName, tParameters[4])   
        
        if ((len(lNucleusNo) == 1 and lNucleusNo[0] == "") 
            or len(lNucleusNo) == 0): 
            lNucleusNoAll = list()
            if tools.CheckExistLocal(sBaseName + ".xyz") == False:        
                print ("\nWARNING! Isotropic sheildings are not exported.\t"
                        + sBaseName + ".xyz" + "\tFile could not be found!")
                return                 
            lNucleusNo = NucleiNoandListXYZ(sBaseName)
            for irow in xrange(0,len(lNucleusNo)):
                lNucleusNoAll.append(lNucleusNo[irow][0])
            lNucleusNo = lNucleusNoAll

        try: 
            iAtomNo = 0
            atom_no = 0
            for i, line in enumerate(lines):
                if 'Atom:' in line:   
                    #atom_no = str(line.split( )[-2])
                    atom_name = str(line.split( )[-1])
                if 'isotropic' in line:    
                    for nucleus_no in range(0, len(tParameters[4])):
                        if atom_name == list(tParameters[4])[nucleus_no]:  
                            atom_no = lNucleusNo[iAtomNo]
                            f2.write(str(atom_no) + '  ' + atom_name)
                            f2.write('  ' + line.split( )[-1] + '\n')	
                    iAtomNo = iAtomNo + 1
                    #if not atom_name in list(tParameters[4]):
                        #try: 
                            #atom_no = int(atom_no) + 1
                        #except:
                            #atom_no = atom_no
                        #f2.write(str(atom_no) + '  ' + atom_name)
                        #f2.write('  ' + '\n')
        except:
            print ("\nWARNING! Isotropic sheildings are not exported.\t"
                        + sBaseName)
        
def chemical_shifts(sheet1, sheet2):
   
    lParameters = scalingparameters(sheet2)

    for i in range(1, sheet2.nrows):
    
        tOriginalParameters, tParameters  = ExtractInputs(sheet2, i)

        if tools.detectformat(tParameters[6]) == "InChI":
            refInChI = tParameters[6]    
            refInChIKey = inchi2input.inchi2inchikey(refInChI)
            sRefBaseName = refInChIKey
        else:
            sRefBaseName = tParameters[6].split(".xyz")[0]    
    
        for irow in range(0, sheet1.nrows):
            
            cell_input = sheet1.cell_value(irow, 0).encode('utf-8')
            inputtype = tools.detectformat(cell_input)
            if inputtype == "InChI":
                InChI, InChIKey = inchi2input.InchIs(sheet1, irow, 0)   
                sBaseName = InChIKey 
            elif inputtype == ".xyz":
                sBaseName = cell_input.split(".xyz")[0]
  
            write_shifts(sBaseName, sRefBaseName, 
                         tOriginalParameters, tParameters)
      
        #print ("\nMethod " + str(i) + " : " 
        #        + tools.method(tOriginalParameters) 
        #        + "\nNMR Chemical Shifts Relative to Reference Molecule of " 
        #        + sRefBaseName 
        #        + " Are Appended to Files With Extension _3D.mol\n")   

        if lParameters[i-1][0] == True:
            
            if tParameters[7] == "YES": 
                
                sBlockName = str("Calculated Shifts" + "\n" 
                                + tools.method(tOriginalParameters))
                lfSlope, lfIntercept, lfR_Squared = (basicoperations.linear_regression(sheet1, 
                                     tOriginalParameters, tParameters, 
                                     sBlockName))    
                    
            else: 

                lScalingShiftInfo = tParameters[7].split(",")
                
                lsSlope = lScalingShiftInfo[:len(tParameters[4])]
                lsIntercept = lScalingShiftInfo[len(tParameters[4]):] 
                lsR_Squared = lScalingShiftInfo[len(tParameters[4]):] 
                
                lfSlope = list()
                lfIntercept = list()
                lfR_Squared = list()
                
                for nucleus_no in xrange(0, len(tParameters[4])):
                   if tools.IsFloat(lsSlope[nucleus_no]) == True:
                       lfSlope.append(float(lsSlope[nucleus_no]))
                   else:
                       lfSlope.append(1.)
                       
                for nucleus_no in xrange(0, len(tParameters[4])):
                   if tools.IsFloat(lsIntercept[nucleus_no]) == True:
                       lfIntercept.append(float(lsIntercept[nucleus_no]))                      
                   else: 
                       lfIntercept.append(0.)

                for nucleus_no in xrange(0, len(tParameters[4])):
                   if tools.IsFloat(lsR_Squared[nucleus_no]) == True:
                      lfR_Squared.append(float(lsR_Squared[nucleus_no]))                      
                   else: 
                      lfR_Squared.append(1.)
                       
            for _irow in xrange(0, sheet1.nrows):
                
                cell_input = sheet1.cell_value(_irow, 0).encode('utf-8')
                inputtype = tools.detectformat(cell_input)
                if inputtype == "InChI":
                    InChI, InChIKey = inchi2input.InchIs(sheet1, _irow, 0)   
                    sBaseName = InChIKey 
                elif inputtype == ".xyz":
                    sBaseName = cell_input.split(".xyz")[0]    
                
                write_scaledshifts(sBaseName, sRefBaseName, 
                                   tOriginalParameters, tParameters, 
                                   lfSlope, lfIntercept, lfR_Squared)       

        if lParameters[i-1][1] == True:

            if tParameters[8] == "YES":   
                
                sBlockName = str("Isotropic Shielding" + "\n" 
                                + tools.method(tOriginalParameters))
                lfSlope, lfIntercept, R_Squared = (
                    basicoperations.linear_regression(sheet1, tOriginalParameters, 
                                                tParameters, sBlockName))

            else: 

                lScalingShiftInfo = tParameters[8].split(",")
                
                lsSlope = lScalingShiftInfo[:len(tParameters[4])]
                lsIntercept = lScalingShiftInfo[len(tParameters[4]):]
                lsR_Squared = lScalingShiftInfo[len(tParameters[4]):] 
                
                lfSlope = list()
                lfIntercept = list()
                lfR_Squared = list()
                
                for ii in xrange(0, len(tParameters[4])):
                   if tools.IsFloat(lsSlope[ii]) == True:
                       lfSlope.append(float(lsSlope[ii]))
                   else:
                       lfSlope.append(1.)
                       
                for ii in xrange(0, len(tParameters[4])):
                   if tools.IsFloat(lsIntercept[ii]) == True:
                       lfIntercept.append(float(lsIntercept[ii]))                      
                   else: 
                       lfIntercept.append(0.)
                for ii in xrange(0, len(tParameters[4])):
                   if tools.IsFloat(lsR_Squared[ii]) == True:
                      lfR_Squared.append(float(lsR_Squared[ii]))                      
                   else: 
                      lfR_Squared.append(1.)
                      
            for _irow in xrange(0, sheet1.nrows):

                cell_input = sheet1.cell_value(_irow, 0).encode('utf-8')
                inputtype = tools.detectformat(cell_input)
                if inputtype == "InChI":
                    InChI, InChIKey = inchi2input.InchIs(sheet1, _irow, 0)   
                    sBaseName = InChIKey 
                elif inputtype == ".xyz":
                    sBaseName = cell_input.split(".xyz")[0]

                write_shifts_regression(sBaseName, tOriginalParameters, 
                                        tParameters, lfSlope, lfIntercept, lfR_Squared)
        
        # print "\nNMR Chemical Shifts by Regression Scaling Are Appended 
        # to Files With Extension _3D.mol\n"   
             
def write_shifts(sBaseName, sRefBaseName, tOriginalParameters, tParameters):

    '''
    This function calculates chemical shifts relative to the 
    reference molecule.  It writes the chemical shifts to the 3D mol 
    files.  It adds the simulation details such as number of nodes, 
    number cores per node and cpu time.  It also adds InChIKey of 
    Reference molecule and the DFT methods used.
    '''

    if tools.CheckExistLocal(sBaseName + "_3D.mol") == False:
        
        print ("\nWARNING! Chemical shifts are not written "
               + "to the file \t" + sBaseName + "_3D.mol" 
               + "\nFile could not be found!")
        
        return
        
    sFileName = (sBaseName + "_"
                    + str(tools.method(tParameters).replace("/","_")) 
                    + ".output")

    if tools.CheckExistLocal(sFileName) == False:
        
        print ("\nWARNING! NMR Chemical Shifts are not "
                  + "calculated for the molecule: "
                  + sBaseName + "\t" + tools.method(tOriginalParameters)
                  + "\t" + sRefBaseName)
        
        return
        
    if tools.CheckExistLocal(sRefBaseName + "_"
                    + str(tools.method(tParameters).replace("/","_")) 
                    + ".output") == False:
        
        print ("\nWARNING! NMR Chemical Shifts are not "
                  + "calculated for the molecule: "
                  + sBaseName + "\t" + tools.method(tOriginalParameters)
                  + "\t" + sRefBaseName + "\tNo such file " + sRefBaseName)
        
        return       
        
    # molecule of ineterest
    if tools.CheckText(sFileName, "Total times") == False:

        print ("\nWARNING! NMR Chemical Shifts are not "
              + "calculated for the molecule: "
              + sBaseName + "\t" + tools.method(tOriginalParameters)
              + "\t" + sRefBaseName)
        return
            
    sFileName = (sRefBaseName + "_"
                + str(tools.method(tParameters).replace("/","_")) 
                + ".output")
    # reference molecule
    if tools.CheckText(sFileName, "Total times") == False:
        
        print ("\nWARNING! NMR Chemical Shifts are not "
              + "calculated for the reference molecule: "
              + sRefBaseName + "\t" + tools.method(tOriginalParameters))
        
        return        
        
    if tools.CheckText(sBaseName + "_3D.mol", 
                       str("Isotropic Shielding" + "\n" 
                           + tools.method(tOriginalParameters))) == False:
        
        print ("\nWARNING! NMR Chemical Shifts are not "
              + "calculated for the molecule: "
              + sBaseName + "\t" + tools.method(tOriginalParameters)
              + "/" + sRefBaseName)
        print ("\nTable is not found in the file: " 
              + sBaseName + "_3D.mol for\n" + "Isotropic Shieldings " 
              + tools.method(tOriginalParameters))
       
        return
         
    if tools.CheckText(sRefBaseName + "_3D.mol", 
                       str("Isotropic Shielding" + "\n" 
                           + tools.method(tOriginalParameters))) == False:
        
        print ("\nWARNING! NMR Chemical Shifts are not "
              + "calculated for the molecule: "
              + sRefBaseName + "\t" + tools.method(tOriginalParameters))
        print ("\nTable is not found in the file: " 
              + sRefBaseName + "_3D.mol for\n" + "Isotropic Shieldings " 
              + tools.method(tOriginalParameters))
        return
            
    dictTable1 = tools.dicTables(sBaseName, 
                                 str("Isotropic Shielding" + "\n" 
                                 + tools.method(tOriginalParameters)))
    dictTables = list() #calculated shifts

    if tools.CheckText(sBaseName + "_3D.mol", 
                       "Calculated Shifts Relative to the Reference of: " 
                       + sRefBaseName + "\n" 
                       + str(tools.method(tOriginalParameters))) == False:
        with open(sBaseName + "_3D.mol", 'a') as f1:
            f1.write("\n\n" 
                    + "Calculated Shifts Relative to the Reference of: " 
                    + sRefBaseName + "\n")
            f1.write(tools.method(tOriginalParameters) + "\n") 
            sFileName = (sBaseName + "_" 
                        + str(tools.method(tParameters).replace("/","_")))
            f1.write("Number of nodes: " 
                    + str(Specs(sFileName)[0]) + "\t" 
                    + "Number of cores per node: " 
                    + str(Specs(sFileName)[1]) + "\t" 
                    + "Total CPU time: " + str(CPUtime(sFileName)) + "\t" 
                    + "Total time: " 
                    + str(Specs(sFileName)[2]) + "s" 
                    + "\n") 
					
            nucleus = list(tParameters[4])
            for i in range(0, len(dictTable1)):
                
                if dictTable1[i][1] in nucleus and len(dictTable1[i]) > 2:
                    for j in range(0, len(tParameters[4])):
                        if dictTable1[i][1] == nucleus[j]:
                            Nuclei_reference_shielding = (
                            basicoperations.average(sRefBaseName, 
                                str("Isotropic Shielding" + "\n" 
                                   + tools.method(tOriginalParameters)), 
                                2, nucleus[j])) 								
                            if sRefBaseName == 'CZDYPVPMEAXLPK-UHFFFAOYSA-N':
                                listNuclei = [dictTable1[i][0], nucleus[j], 
                                str(round(Nuclei_reference_shielding
                                - float(dictTable1[i][2]),4))]
                            else:					
                                Nuclei_reference_expshift = (
                                basicoperations.average(sRefBaseName, 
                                    "Experimental Shifts", 2, nucleus[j])) 
                                listNuclei = [dictTable1[i][0], nucleus[j], 
                                str(round(Nuclei_reference_shielding 
                                - float(dictTable1[i][2]) 
                                - Nuclei_reference_expshift,4))]
                            dictTables.append(listNuclei)               
                else: 
                    listOthers = [dictTable1[i][0], dictTable1[i][1]]
                    dictTables.append(listOthers)   
            
            for listLine in dictTables:
                for sWord in listLine:
                    f1.write(sWord + "\t")
                f1.write("\n") 
                
def write_scaledshifts(sBaseName, sRefBaseName, tOriginalParameters, 
                       tParameters, lSlope, lIntercept, lR_Squared):
    
    '''
    This function calculates chemical shifts by linear regression.
    It writes the chemical shifts to the 3D mol files.
    It writes the simulation details such as the number of nodes, 
    number of cores per node and cpu time.
    It also writes Slope, Intercept and the DFT methods used 
    "Computed Chemical Shifts vs. Experimental Shifts"
    '''
    
    sFileName = (sRefBaseName + "_"
                        + str(tools.method(tParameters).replace("/","_")) 
                        + ".output")  
    
    if tools.CheckExistLocal(sFileName) == False:
      
      print ("\nWARNING! NMR Chemical Shifts are not calculated "
               + "for the molecule: " + sBaseName + "\tMethod: " 
               + tools.method(tOriginalParameters) + "/" + sRefBaseName)
      
      return
         
    if tools.CheckText(sFileName, "Total times") == False:
        
        print ("\nWARNING! NMR Chemical Shifts are not "
                  + "calculated for the molecule: " + sRefBaseName + "\t" 
                  + tools.method(tOriginalParameters))
        
        return
        
    sFileName = (sBaseName + "_"
                    + str(tools.method(tParameters).replace("/","_")) 
                    + ".output")
    
    if tools.CheckExistLocal(sFileName) == False:
            
        print ("\nWARNING! NMR Chemical Shifts are not calculated "
               + "for the molecule: " + sRefBaseName + "\tMethod: " 
               + tools.method(tOriginalParameters))
      
        
        return    
    
    if tools.CheckText(sFileName, "Total times") == False:
        
        print ("\nWARNING! NMR Chemical Shifts are not "
                  + "calculated for the molecule: " + sBaseName 
                  + "\tMethod: " + tools.method(tOriginalParameters)
                  + "/" + sRefBaseName)
            
        return
             
    dictTable1 = tools.dicTables(sBaseName, 
                                 str("Calculated Shifts Relative to the " 
                                 + "Reference of: " + sRefBaseName + "\n" 
                                 + tools.method(tOriginalParameters)))
    dictTables = list() #calculated shifts
    
    if (tools.CheckExistLocal(sBaseName + "_3D.mol") == True and 
        tools.CheckText(sBaseName + "_3D.mol", 
            "Scaled Shifts Relative to the Reference of: " 
                      + sRefBaseName + "\n" 
                      + str(tools.method(tOriginalParameters))) == False):
        
        with open(sBaseName + "_3D.mol", 'a') as f1:
            f1.write("\n\n" 
                    + "Scaled Shifts Relative to the Reference of: " 
                    + sRefBaseName 
                    + "\n"
                    + tools.method(tOriginalParameters) 
                    + "\n") 
            sFileName = (sBaseName + "_" 
                        + str(tools.method(tParameters).replace("/","_")))
            f1.write("Number of nodes: " 
                    + str(Specs(sFileName)[0]) + "\t" 
                    + "Number of cores per node: " 
                    + str(Specs(sFileName)[1]) + "\t" 
                    + "Total CPU time: " + str(CPUtime(sFileName)) + "\t" 
                    + "Total time: " 
                    + str(Specs(sFileName)[2]) + "s" 
                    + "\n")
            f1.write("Slope: " + str(lSlope) + "\t" + "Intercept: " 
                    + str(lIntercept) + "\t" + "R_squared: " + str(lR_Squared) + "\n")
            
            for i in range(0, len(dictTable1)):
                
                nucleus = list(tParameters[4])

                if dictTable1[i][1] in nucleus and len(dictTable1[i])>2:
                    for nucleus_no in range(0, len(tParameters[4])):
                        if dictTable1[i][1] == nucleus[nucleus_no]: 
                            listNucleus = [dictTable1[i][0], 
                                          nucleus[nucleus_no], 
                                          str(round((float(dictTable1[i][2]) - 
                                                      lIntercept[nucleus_no])
                                                      /lSlope[nucleus_no],4))]
                            dictTables.append(listNucleus)
                else:
                    listOthers = [dictTable1[i][0], dictTable1[i][1]]
                    dictTables.append(listOthers) 

            #f1.write(dictTables)
            for listLine in dictTables:
                for sWord in listLine:
                    f1.write(sWord + "\t")
                f1.write("\n")
            #print dictTables
                        
def write_shifts_regression(sBaseName, tOriginalParameters, tParameters, 
                            SlopeList, InterceptList, R_SquaredList):

    "Isotropic shieldings vs. Experimental Shifts"
    
    sFileName = (sBaseName + "_"
                    + str(tools.method(tParameters).replace("/","_")) 
                    + ".output")
    
    if tools.CheckExistLocal(sFileName) == False:
        
        print ("\nWARNING! NMR Chemical Shifts are not calculated for the molecule: "
               + sBaseName + "\t" + tools.method(tOriginalParameters))    
        return  
        
    if tools.CheckText(sFileName, "Total times") == False:
        
        print ("\nWARNING! NMR Chemical Shifts are not calculated for the molecule: "
                  + sBaseName + "\t" + tools.method(tOriginalParameters))
        return
        
    dictTable1 = tools.dicTables(sBaseName, str("Isotropic Shielding" + "\n" 
                                + tools.method(tOriginalParameters)))
    dictTables = list() #calculated shifts

    if tools.CheckText(sBaseName + "_3D.mol", 
                       "Predicted Shifts by Linear Regression" + "\n" 
                       + str(tools.method(tOriginalParameters))) == False:
        with open(sBaseName + "_3D.mol", 'a') as f1:
            f1.write("\n\nPredicted Shifts by Linear Regression\n")
            f1.write(tools.method(tOriginalParameters) + "\n") 
            sFileName = (sBaseName + "_" 
                        + str(tools.method(tParameters).replace("/","_")))
            f1.write("Number of nodes: " 
                    + str(Specs(sFileName)[0]) + "\t" 
                    + "Number of cores per node: " 
                    + str(Specs(sFileName)[1]) + "\t" 
                    + "Total CPU time: " + str(CPUtime(sFileName)) + "\t" 
                    + "Total time: " 
                    + str(Specs(sFileName)[2]) + "s" 
                    + "\n")
            f1.write("Slope: " + str(SlopeList) + "\t" + "Intercept: " 
                    + str(InterceptList) + "\t" + "R_Squared: " + str(R_SquaredList) + "\n")
            
            for i in range(0, len(dictTable1)):
                
                nucleus = list(tParameters[4])
                
                if dictTable1[i][1] in nucleus and len(dictTable1[i])>2:
                    for nucleus_no in range(0, len(tParameters[4])):
                        if dictTable1[i][1] == nucleus[nucleus_no]:  
                            listNucleus = [dictTable1[i][0], 
                                          nucleus[nucleus_no], 
                                          str(round((float(dictTable1[i][2]) - 
                                              InterceptList[nucleus_no])
                                              /SlopeList[nucleus_no],4))]
                            dictTables.append(listNucleus)		              
                else: 
                    listOthers = [dictTable1[i][0], dictTable1[i][1]]
                    dictTables.append(listOthers)      

            #f1.write(dictTables)
            for listLine in dictTables:
                for sWord in listLine:
                    f1.write(sWord + "\t")
                f1.write("\n")
            #print dictTables   
               
def scalingparameters(sheet2):
    
    lParameters = list()
    
    for irow in range(1, sheet2.nrows):
    
        tOriginalParameters, tParameters  = ExtractInputs(sheet2, irow)

        if tParameters[7] == "YES": bDoCal1 = True
        else: 
            bDoCal1 = True
            for l in tParameters[7].split(","):
                if tools.IsFloat(l) == False:
                        bDoCal1 = False
        
        if tParameters[8] == "YES": bDoCal2 = True
        else: 
            bDoCal2 = True
            for l in tParameters[8].split(","):
                if tools.IsFloat(l) == False:
                        bDoCal2 = False
                                 
        lParameters.append([bDoCal1, bDoCal2])
    
    return lParameters
    
def CPUtime(sFileName):
    
    '''This function extracts total cpu times from a given output file.'''
    
    sCputime = ''
    for name in glob.glob(sFileName + '.output*'):
        with open(name) as f1:
            lines = f1.readlines() 
            for i, line in enumerate(lines):
                if 'Total times  cpu:' in line:   			
                    sCputime = str(line.split( )[3])	
                
    return sCputime
 
def extractNODE():

    try:
        with open('finalsetrun.bash', 'r') as file:
            lines = file.readlines()
        for i, line in enumerate(lines):
                if 'srun -N ' in line:   			
                    return str(line.split("srun -N ")[-1].split(" ")[0])	     
    except:
        return str(0)

def extractPROCESSOR():

    try:
        with open('finalsetrun.bash', 'r') as file:
            lines = file.readlines()
        for i, line in enumerate(lines):
                if 'srun -N ' in line:   			
                    return str(line.split("--exclusive --mpi=pmi2 -n ")[-1].split(" ")[0])	     
    except:
        return str(0)
        
def Specs(sFileName):
    
    '''This function gives the simulation specs'''    
    
    fCputime = (float(CPUtime(sFileName).split("s")[0]) 
                    if "s" in CPUtime(sFileName) 
                    else float(CPUtime(sFileName)))

    sNODES = extractNODE()
    sPROCESSOR = extractPROCESSOR()
    
    return sNODES, sPROCESSOR, fCputime   
    
def shielding_index(sBaseName, sNuclei):
    
    '''
    This function takes the list of nuclei which their isotropic 
    shiledings will be calculated.  It gives the total number of 
    nuclei and indexes of each nucleus.
    '''
    
    lNuclei = ReadXYZFile(sBaseName + ".xyz")[1]
    lNucleusNo = list()    
    
    for no in range(0, len(list(sNuclei))):
        for i in range(0, len(lNuclei)):
            if list(sNuclei)[no] == lNuclei[i]:
                lNucleusNo.append(i)    

    sNucleusNo = ""
    for i in range(0, len(lNucleusNo)):
        sNucleusNo = sNucleusNo + str(lNucleusNo[i] + 1) + " "

    return str(len(lNucleusNo)) + " " + sNucleusNo

def read_shielding_index(sFileName):
    
    for name in glob.glob(sFileName):
        with open(name) as f1:
            lines = f1.readlines() 
            
    lsIndecies = list()
         
    for i, line in enumerate(lines):
        if "SHIELDING" in line:   
            line1 = line.split("SHIELDING")[-1]
            lsIndecies = line1.split(" ")
    
    lsIndecies = [x for x in lsIndecies if x != ""]
    lsIndecies = [x.replace("\n","") if "\n" in x else x for x in lsIndecies]
    
    try:
        lsIndecies[0] = int(lsIndecies[0])
    except:
        lsIndecies[0] = lsIndecies[0]
        
    if type(lsIndecies[0]) == int and lsIndecies[0] == len(lsIndecies)-1:
        del lsIndecies[0]
            
    return lsIndecies

def NucleiNoandListXYZ(sBaseName):
    
    '''
    This function lists the number of each atom an its symbol 
    by reading the xyz file
    '''
    lsXYZAtomMolecule = ReadXYZFile(sBaseName + ".xyz")[1]  
    
    lNuclei = list()
    
    iCurrAtom = 1
    for irow in range(0, len(lsXYZAtomMolecule)):
        ListNucleus = iCurrAtom, lsXYZAtomMolecule[irow]
        lNuclei.append(ListNucleus)
        iCurrAtom = iCurrAtom + 1
        
    return lNuclei

def NucleusNoListXYZ(sBaseName, sNuclei):

    '''This function lists the number of each atom of a molecule'''
    
    if tools.CheckExistLocal(sBaseName + ".xyz") == True:
        
        lsXYZAtomMolecule = ReadXYZFile(sBaseName + ".xyz")[1] 
        lNucleusNo = list()
    
        for i in range(0, len(lsXYZAtomMolecule)):
            for no in range(0, len(list(sNuclei))):           
                if lsXYZAtomMolecule[i] == list(sNuclei)[no]:
                    lNucleusNo.append(str(i+1))
    
        return lNucleusNo
        
    else: 
        
        return 0
    
def ExtractInputs(sheet2, i):
    
    '''
    This method extracts the DFT methods listed in the Excel file and returns 
    them in a parameter list with the order below:
    1: Functional for chemical shift calculation
    2: Basis set for chemical shift calculation
    3: Functional for geometry optimization
    4: Basis set for geometry optimization 
    5: Solvent 
    6: Nuclei list
    7: InChI code of Reference molecule
    8: Scaling option by using calculated shifts
    9: Scaling option by using calculated shieldings
    '''
    
    sNmrOriginalFunctional = sheet2.cell_value(i,0).encode('utf-8')
    if "*" in sNmrOriginalFunctional:
        sNmrFunctional = sNmrOriginalFunctional.replace("*","STAR")
    elif " " in sNmrOriginalFunctional:
        sNmrFunctional = sNmrOriginalFunctional.replace(" ","_")
    else: 
        sNmrFunctional = sheet2.cell_value(i,0).encode('utf-8')

    sNmrOriginalBasis = sheet2.cell_value(i,1).encode('utf-8')
    if "*" in sNmrOriginalBasis:
        sNmrBasis = sNmrOriginalBasis.replace("*","STAR")
    elif " " in sNmrOriginalBasis:
        sNmrBasis = sNmrOriginalBasis.replace(" ","_")
    elif "+" in sNmrOriginalBasis:
        sNmrBasis = sNmrOriginalBasis.replace(" ","PLUS")
    elif "(" in sNmrOriginalBasis:
        sNmrBasis = sNmrOriginalBasis.replace(" ","")
    elif ")" in sNmrOriginalBasis:
        sNmrBasis = sNmrOriginalBasis.replace(" ","")
    else: 
        sNmrBasis = sheet2.cell_value(i,1).encode('utf-8')

    sGeometryOriginalFunctional = sheet2.cell_value(i,2).encode('utf-8')
    if "*" in sGeometryOriginalFunctional:
        sGeometryFunctional = sNmrOriginalFunctional.replace("*","STAR")
    elif " " in sGeometryOriginalFunctional:
        sGeometryFunctional = sGeometryOriginalFunctional.replace(" ","_")
    else: 
        sGeometryFunctional = sheet2.cell_value(i,2).encode('utf-8')
        
    sGeometryOriginalBasis = sheet2.cell_value(i,3).encode('utf-8')
    if "*" in sGeometryOriginalBasis:
        sGeometryBasis = sGeometryOriginalBasis.replace("*","STAR")
    elif " " in sGeometryOriginalBasis:
        sGeometryBasis = sGeometryOriginalBasis.replace(" ","_")
    elif "+" in sGeometryOriginalBasis:
        sGeometryBasis = sGeometryOriginalBasis.replace(" ","PLUS")
    elif "(" in sGeometryOriginalBasis:
        sGeometryBasis = sGeometryOriginalBasis.replace(" ","")
    elif ")" in sGeometryOriginalBasis:
        sGeometryBasis = sGeometryOriginalBasis.replace(" ","")
    else: 
        sGeometryBasis = sheet2.cell_value(i,3).encode('utf-8')
        
    sNuclei = sheet2.cell_value(i, 4).encode('utf-8').upper()
    
    if (sheet2.cell_value(i, 5).encode('utf-8')).upper() == "GAS":
        sSolvent = sheet2.cell_value(i, 5).encode('utf-8').upper()
    else:
        sSolvent = sheet2.cell_value(i, 5).encode('utf-8')

    sRefBaseName = sheet2.cell_value(i, 6).encode('utf-8')

    sScalingShift = sheet2.cell_value(i, 7).encode('utf-8').upper()
    sScalingShielding = sheet2.cell_value(i, 8).encode('utf-8').upper()

    tParameters = (sNmrFunctional, sNmrBasis, sGeometryFunctional, 
                   sGeometryBasis, sNuclei, sSolvent, sRefBaseName, 
                   sScalingShift, sScalingShielding)
    tOriginalParameters = (sNmrOriginalFunctional, sNmrOriginalBasis, 
                           sGeometryOriginalFunctional, sGeometryOriginalBasis, 
                           sNuclei, sSolvent, sRefBaseName, sScalingShift, 
                           sScalingShielding)

    return tOriginalParameters, tParameters #type: tuple

def ArrangeInputs(sLocalDir, sServerDir, sFileName1, sFileName2):
    
    sLocalDir = sLocalDir.replace("\\","/")
    sLocalDir = sLocalDir.replace("\t","/t")
    sLocalDir = (sLocalDir + "/" if sLocalDir.split("/")[-1] != "" 
                    else sLocalDir)
    sServerDir = sServerDir.replace("\\","/")
    sServerDir = sServerDir.replace("\t","/t")
    sServerDir = (sServerDir + "/" if (sServerDir.split("/")[-1] != "") 
                    else sServerDir)
    sServerDir = ("/" + sServerDir if (sServerDir.split("/")[0] != "") 
                    else sServerDir)
    
    # takes the excel sheet # Molecule file
    workbook1 = xlrd.open_workbook(sFileName1) 
    # Method file
    workbook2 = xlrd.open_workbook(sFileName2) 
    # takes 1st tab # Molecule set
    sheet1 = workbook1.sheet_by_index(0) 
    # Method file
    sheet2 = workbook2.sheet_by_index(0) 

    return sLocalDir, sServerDir, sheet1, sheet2

def ReadXYZFile(sFileName):
    '''
    This method reads a given xyz file. Filename includes the extension of .xyz. 
    Ex: UHOVQNZJYSORNB-UHFFFAOYSA-N.xyz  It skips the first two lines. 
    First line is the total number of atoms in the file. Second line 
    indicates the program which the file was generated by.  
    
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
    lsXYZFile = list()

    with open(sFileName, 'r') as f:    
            lines = f.readlines() 
            for i, line in enumerate(lines):
                # It skips the first two lines 
                # First line: the number of atoms in the file
                # First line: + generated by VMD 
                if i >= 2 and len(line.split())>0: 
                    lsXYZFile.append(line.split()) 

    lsXYZCoorMolecule = [(lsXYZFile[row][1:]) 
                        for row in xrange(0,len(lsXYZFile))]
    lsXYZAtomMolecule = [(lsXYZFile[row][0]) 
                        for row in xrange(0,len(lsXYZFile))]
    lfXYZCoorMolecule = list()
    for row in xrange(0,len(lsXYZFile)):
        for column in xrange(0,len(lsXYZCoorMolecule[row])):
            try: 
                lfXYZCoorMolecule.append(float(lsXYZCoorMolecule[row][column])) 
            except:
                lfXYZCoorMolecule.append(None) 
    lfXYZCoorMolecule = numpy.reshape(lfXYZCoorMolecule, 
                                     (len(lsXYZCoorMolecule), -1), 
                                     order='F')
     
    # return type: list                 
    return lsXYZFile, lsXYZAtomMolecule, lsXYZCoorMolecule, lfXYZCoorMolecule

  