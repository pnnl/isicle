# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 11:24:24 2016

@author: yesi172
"""

import math
import numpy

import inchi2input
import tools    
    
def average(sFileName, sBlockName, j, sNuclei):
    
    '''
    Function to calculate the average value of isotropic shieldings 
    or chemical shifts of a molecule for a given nucleus.
    sFileName is the file which consists of tables. Each table 
    is represented by a sBlockName (i.e. Experimental Shifts, 
    Isotropic Shieldings, ..). Each block consists of a series of nuclei 
    and their corresponding data.
    '''
    
    dictTable = tools.dicTables(sFileName, sBlockName) 
    
    Nucleus = 0 
    count_Nucleus = 0

    for i in range(0,len(dictTable)):	
        if dictTable[i][1] == sNuclei and len(dictTable[i])>2:
            Nucleus = Nucleus + float(dictTable[i][j])
            count_Nucleus = count_Nucleus + 1
            
    if count_Nucleus != 0:    
        return Nucleus / count_Nucleus
    else: 
        return ""                   

def average_col(sheet, iTopRow, iBottomRow, iColumn):
    
    '''
    This function calculates the average of a column between the 
    top and bottom row of a given Excel sheet
    ''' 
    
    cell_value = 0.
    count = 0.

    for irow in xrange(iTopRow, iBottomRow+1):
        try:
            cell_value = cell_value + float(sheet.cell_value(irow, iColumn))
            count = count + 1
        except:
            continue

    if count != 0.:    
        return cell_value / count
    else: 
        return ""
 
def mae(sBaseName, sBlock, sNuclei):
    
    '''
    This function calculates the mean absolute error of the given nuclei 
    or nucleus in a given molecule.  sBlock represents the 
    table on 3D mol files.
    '''  
    
    try:
        dictTable1 = tools.dicTables(sBaseName, sBlock)
    except:
        return
    
    try:
        dictTable2 = tools.dicTables(sBaseName, "Experimental Shifts")
    except:
        return

    Nucleus = 0. 
    countNucleus = 0.

    for i in range(0, len(dictTable1)):
        j = int(dictTable1[i][0])-1
        try:
            if (dictTable1[i][1] == sNuclei and len(dictTable1[i]) > 2 and 
                    len(dictTable2[j]) > 2 and dictTable2[j][1] == sNuclei):
                Nucleus = (Nucleus + abs(float(dictTable1[i][2]) 
                                        - float(dictTable2[j][2])))
                countNucleus = countNucleus + 1
        except: 
            continue

    if countNucleus != 0:     
        return Nucleus/countNucleus
    else: 
        return ""

def mse(sBaseName, sBlock, sNuclei):
    
    '''
    This function calculates mean square error of the nuclei or nucleus 
    in a given molecule.  sBlock represents the table that 
    will be used.
    ''' 
        
    dictTable1 = tools.dicTables(sBaseName, sBlock)
    dictTable2 = tools.dicTables(sBaseName, "Experimental Shifts")
    
    mseNucleus = 0 
    countNucleus = 0
    
    for i in range(0,len(dictTable1)):
        j = int(dictTable1[i][0])-1
        try:
            if (dictTable1[i][1] == sNuclei and len(dictTable1[i]) > 2 and 
                    len(dictTable2[j]) > 2 and dictTable2[j][1] == sNuclei):
                mseNucleus = mseNucleus + math.pow(abs(float(dictTable1[i][2]) 
                                        - float(dictTable2[j][2])),2)    
                countNucleus = countNucleus + 1
        except:
            continue
        
    if countNucleus != 0:     
        return float(mseNucleus/countNucleus)
    else: 
        return ""    
    
def rmse(sBaseName, sBlock, sNuclei):
    
    '''
    This function calculates root mean square error of the nuclei or 
    nucleus in a given molecule.  sBlock represents the table 
    that will be used.
    '''  
    try:     
        return math.sqrt(mse(sBaseName, sBlock, sNuclei))
    except: 
        return ""

def msie(sBaseName, sBlock, sNuclei):   
    
    '''    
    This function calculates mean signed error of the nuclei or nucleus 
    in a given molecule as InChIKey. sBlock represents the table that 
    will be used.
    '''
    
    dictTable1 = tools.dicTables(sBaseName, sBlock)
    dictTable2 = tools.dicTables(sBaseName, "Experimental Shifts")
    
    Nucleus = 0 
    countNucleus = 0
    
    for i in range(0,len(dictTable1)):
        j = int(dictTable1[i][0])-1
        try:
            if (dictTable1[i][1] == sNuclei and len(dictTable1[i]) > 2 and 
                    len(dictTable2[j]) > 2 and dictTable2[j][1] == sNuclei):
                Nucleus = (Nucleus + float(dictTable1[i][2]) 
                                    - float(dictTable2[j][2]))
                countNucleus = countNucleus + 1                     
        except:
            continue
        
    if countNucleus != 0:    
        return Nucleus/countNucleus
    else: 
        return ""              

def maxae(sBaseName, sBlock, sNuclei):
    
    '''    
    This function calculates the maximum absolute error of the nuclei 
    or nucleus in a given molecule. sBlock represents the 
    table that will be used.
    ''' 
    
    dictTable1 = tools.dicTables(sBaseName, sBlock)
    dictTable2 = tools.dicTables(sBaseName, "Experimental Shifts")
    
    listNucleus = list()
    countNucleus = 0
    
    for i in range(0,len(dictTable1)):	
        j = int(dictTable1[i][0])-1
        try:
            if (dictTable1[i][1] == sNuclei and len(dictTable1[i]) > 2 and 
                    len(dictTable2[j]) > 2 and dictTable2[j][1] == sNuclei):
                (listNucleus.append(abs(float(dictTable1[i][2]) 
                                        - float(dictTable2[j][2]))))
                countNucleus = countNucleus + 1     
        except:
            continue
        
    if countNucleus != 0:    
        return max(listNucleus)
    else: 
        return ""

def exp_compt_shifts(sBaseName, sBlock, sNuclei, lBondedAtomNumber):
    
    '''    
    This function calculates the absolute error of the nuclei 
    or nucleus in a given molecule. sBlock represents the 
    table that will be used.
    ''' 
    
    dictTable1 = tools.dicTables(sBaseName, sBlock)
    dictTable2 = tools.dicTables(sBaseName, "Experimental Shifts")
    
    listNucleusExp = list()
    listNucleusComp = list()
    listAtomName = list()
    listAtomNumber = list()
    
    for i in range(0,len(dictTable1)):	
        for k in xrange(0,len(lBondedAtomNumber)):
            j = int(dictTable1[i][0])-1
            try:
                if (dictTable1[i][1] == sNuclei and len(dictTable1[i]) > 2 and 
                        len(dictTable2[j]) > 2 and dictTable2[j][1] == sNuclei and
                        lBondedAtomNumber[k] == dictTable1[i][0] and
                        lBondedAtomNumber[k] == dictTable2[j][0]):
                    listNucleusComp.append(float(dictTable1[i][2]))
                    listNucleusExp.append(float(dictTable2[j][2]))
                    listAtomName.append(dictTable1[i][1])
                    listAtomNumber.append(dictTable1[i][0])                
            except:
                continue
        
    listErrors = zip(listAtomNumber, listAtomName, listNucleusExp, listNucleusComp)
    
    return listErrors

def error_sorted(sBaseName, sBlock, sNuclei, order):
    
    '''    
    This function calculates the absolute error of the nuclei 
    or nucleus in a given molecule. sBlock represents the 
    table that will be used.
    ''' 
    
    dictTable1 = tools.dicTables(sBaseName, sBlock)
    dictTable2 = tools.dicTables(sBaseName, "Experimental Shifts")
    
    listNucleus = list()
    listAtomName = list()
    listAtomNumber = list()
    
    for i in range(0,len(dictTable1)):	
        j = int(dictTable1[i][0])-1
        try:
            if (dictTable1[i][1] == sNuclei and len(dictTable1[i]) > 2 and 
                    len(dictTable2[j]) > 2 and dictTable2[j][1] == sNuclei):
                (listNucleus.append(abs(float(dictTable1[i][2]) 
                                        - float(dictTable2[j][2]))))
                listAtomName.append(dictTable1[i][1])
                listAtomNumber.append(dictTable1[i][0])                
        except:
            continue
        
    listErrors = zip(listAtomNumber, listAtomName, listNucleus)
    
    if order == 'ascending':
        return sorted(listErrors, key=lambda tup: tup[2]) 
    elif order == 'decending':
        return sorted(listErrors, key=lambda tup: tup[2], reverse=True)
    else:
        return listErrors
        
def std(sheet, iTopRow, iBottomRow, iColumn):

    '''
    This function gives the standard deviation of a column between the 
    given top and bottom rows of an Excel sheet
    ''' 
    
    lCells = list()
    count = 0

    for irow in xrange(iTopRow, iBottomRow+1):
        try:
            lCells.append(float(sheet.cell_value(irow, iColumn)))
            count = count + 1
        except:
            continue

    if count != 0:    
        return numpy.std(lCells)
    else: 
        return ""

def _polyfit(x, y, degree):

    #slope, intercept
    coeffs = numpy.polyfit(x, y, degree)

     # Polynomial Coefficients
    polynomial_coefficients = coeffs.tolist()
    slope = polynomial_coefficients[0]
    intercept = polynomial_coefficients[1]
    
    # r-squared
    p = numpy.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = numpy.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = numpy.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = numpy.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    rvalue = ssreg / sstot

    return slope, intercept, rvalue
    
def linear_regression(sheet, tOriginalParameters, tParameters, sBlockName):
    
    '''
    This function does linear regression scaling to the predicted 
    chemical shifts.  It doesscaling to either isotropic shieldings vs. 
    experimental shifts or calculated shifts to experimental shifts.
    ''' 
    
    SlopeList = list()
    InterceptList = list()  
    R_SquaredList = list()    
    
    for nucleus_no in range(0, len(tParameters[4])):

        listNucleusx = list()
        listNucleusy = list()
              
        for row in range(0, sheet.nrows):
            cell_input = sheet.cell_value(row, 0).encode('utf-8')
            inputtype = tools.detectformat(cell_input)
            if inputtype == "InChI":
                InChI, InChIKey = inchi2input.InchIs(sheet, row, 0)   
                sBaseName = InChIKey 
            elif inputtype == ".xyz":
                sBaseName = cell_input.split(".xyz")[0]
                 
            dictExperimentalShift = (tools.dicTables(sBaseName, 
                                                     "Experimental Shifts"))
            dictTable = tools.dicTables(sBaseName, sBlockName) 

            if dictExperimentalShift is not None and dictTable is not None:
                for i in xrange(0, len(dictExperimentalShift)):   
                    for j in xrange(0, len(dictTable)):
                        if (dictExperimentalShift[i][1] == 
                            list(tParameters[4])[nucleus_no] and 
                            len(dictExperimentalShift[i]) > 2 and
                            len(dictTable[j]) > 2 and 
                            dictExperimentalShift[i][0] == dictTable[j][0]):   
                            # Nucleus #x axis
                            listNucleusx.append(float(dictExperimentalShift[i][2]))
                            # Nucleus #y axis 
                            listNucleusy.append(float(dictTable[j][2]))
                            print sBaseName, i, j,  dictExperimentalShift[i][0], dictTable[j][0], float(dictExperimentalShift[i][2]), float(dictTable[j][2]), float(dictExperimentalShift[i][2]) - float(dictTable[j][2]) 
        
        try: 
            # slope, intercept, r_value, p_value, std_err = 
            # stats.linregress(listNucleusx, listNucleusy) #Nucleus
            slope, intercept, rvalue = _polyfit(listNucleusx, listNucleusy, 1) 
            
        except:
            slope = 1. 
            intercept = 0.
            rvalue = 1.
           
        SlopeList.append(slope)   
        InterceptList.append(intercept)
        R_SquaredList.append(rvalue)
        
    return SlopeList, InterceptList, R_SquaredList

def walltime_str2int(sWalltime):
    
    '''This function converts the time string to integer type.'''
    
    return (int(sWalltime.split(":")[0])*3600 
            + int(sWalltime.split(":")[1])*60 
            + int(sWalltime.split(":")[2]))
            