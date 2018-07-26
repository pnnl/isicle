# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 11:24:24 2016

@author: yesi172
"""

import shutil
import xlrd
import xlsxwriter
import math

import inchi2input
import basicoperations
import tools
import exporters
                  
def summary_file(SHEET1, SHEET2):

    '''It prepares the Summary File (Excel)'''    
     
    lParameters = exporters.scalingparameters(SHEET2)

    workbook = xlsxwriter.Workbook("Summaries.xlsx")
    # Take the excel sheet 
    sheet = workbook.add_worksheet()   
    
    sheet.write(0, 0, "Methods") # row, column, item    
    
    iSummaryErrors = 3
    
    sLongestNuclei = ""
    for irow in xrange(1, SHEET2.nrows):  
        tOriginalParameters, tParameters = exporters.ExtractInputs(SHEET2, irow)
        sNuclei = tParameters[4]
        for inucleus in range(0, len(sNuclei)):    
            if not sNuclei[inucleus] in sLongestNuclei:
                sLongestNuclei = sLongestNuclei + sNuclei[inucleus] 

    sLongestNucleiShift = ""
    for irow in xrange(1, SHEET2.nrows):  
        tOriginalParameters, tParameters = exporters.ExtractInputs(SHEET2, irow)
        sNuclei = tParameters[4]
        if lParameters[irow-1][0] == True: #if tParameters[7] == "YES":
            for inucleus in range(0, len(sNuclei)):    
                if not sNuclei[inucleus] in sLongestNucleiShift:
                    sLongestNucleiShift = sLongestNucleiShift + sNuclei[inucleus]    

    sLongestNucleiShielding = ""
    for irow in xrange(1, SHEET2.nrows):  
        tOriginalParameters, tParameters = exporters.ExtractInputs(SHEET2, irow)
        sNuclei = tParameters[4]
        if lParameters[irow-1][1] == True: #if tParameters[8] == "YES":
            for inucleus in range(0, len(sNuclei)):    
                if not sNuclei[inucleus] in sLongestNucleiShielding:
                    sLongestNucleiShielding = sLongestNucleiShielding + sNuclei[inucleus]

    ShiftBar = 0
    ShieldingBar = 0    
    for irow in xrange(1, SHEET2.nrows):         
        if lParameters[irow-1][0] == True:
             ShiftBar = ShiftBar + 1          
        if lParameters[irow-1][1] == True:
             ShieldingBar = ShieldingBar + 1


    iRightColumn = len(sLongestNuclei)*iSummaryErrors    
    sheet.merge_range(1, 1, 1, len(sLongestNuclei)*iSummaryErrors, 
                      "Errors of Shifts Calculated Relative to"
                      + " the Reference of " 
                      + str(tParameters[6]))     
    for i in xrange(0, len(sLongestNuclei)):
        sheet.write(2, 1+i*iSummaryErrors, 
                    "Mean Absolute Error for " + list(sLongestNuclei)[i])
        sheet.write(2, 2+i*iSummaryErrors, 
                    "Root Mean Square Error for " + list(sLongestNuclei)[i])
        sheet.write(2, 3+i*iSummaryErrors, 
                    "Maximum Absolute Error for " + list(sLongestNuclei)[i])    
    
    if ShiftBar >= 1: 
        iRightColumn = (len(sLongestNuclei)*iSummaryErrors 
                       + len(sLongestNucleiShift)*iSummaryErrors)

        if len(sLongestNucleiShielding) > 0:
            sheet.merge_range(1, len(sLongestNuclei)*iSummaryErrors + 1, 
                              1, len(sLongestNuclei)*iSummaryErrors + 1
                                  + len(sLongestNucleiShift)*iSummaryErrors - 1, 
                             "Errors of Scaled Shifts Relative to the"
                                 + " Reference of " + str(tParameters[6])) 
        else:
            sheet.merge_range(1, len(sLongestNuclei)*iSummaryErrors + 1, 
                              1, len(sLongestNuclei)*iSummaryErrors + 1 
                                  + len(sLongestNucleiShift)*iSummaryErrors, 
                             "Errors of Scaled Shifts Relative to the"
                                 + " Reference of " + str(tParameters[6])) 
        for i in xrange(0, len(sLongestNucleiShift)):
            sheet.write(2, 
                        1+iSummaryErrors*len(sLongestNuclei)+i*iSummaryErrors, 
                        "Mean Absolute Error for " + list(sLongestNuclei)[i])
            sheet.write(2, 
                       2+iSummaryErrors*len(sLongestNuclei)+i*iSummaryErrors, 
                       "Root Mean Square Error for " + list(sLongestNuclei)[i])
            sheet.write(2, 
                       3+iSummaryErrors*len(sLongestNuclei)+i*iSummaryErrors, 
                       "Maximum Absolute Error for " + list(sLongestNuclei)[i])     
    
    if ShiftBar == 0 and ShieldingBar >= 1: 
        iRightColumn = (len(sLongestNuclei)*iSummaryErrors 
                       + len(sLongestNucleiShielding)*iSummaryErrors 
                       if len(sLongestNuclei)*iSummaryErrors 
                           + len(sLongestNucleiShielding)*iSummaryErrors > 
                           iRightColumn 
                       else iRightColumn)
        if len(sLongestNucleiShielding) > 0: 
            sheet.merge_range(1, len(sLongestNuclei)*iSummaryErrors + 1, 
                              1, len(sLongestNuclei)*iSummaryErrors + 1
                                 + len(sLongestNucleiShielding)*iSummaryErrors - 1, 
                              "Errors of Scaled Shifts Shielding")
        else:    
            sheet.merge_range(1, len(sLongestNuclei)*iSummaryErrors + 1, 
                              1, len(sLongestNuclei)*iSummaryErrors + 1
                                  + len(sLongestNucleiShielding)*iSummaryErrors, 
                              "Errors of Scaled Shifts Shielding")
        for i in xrange(0,len(sLongestNucleiShielding)):
            sheet.write(2, 1 + iSummaryErrors*len(sLongestNuclei) 
                           + i*iSummaryErrors, 
                        "Mean Absolute Error for " + list(sLongestNuclei)[i])
            sheet.write(2, 2 + iSummaryErrors*len(sLongestNuclei) 
                           + i*iSummaryErrors, 
                       "Root Mean Square Error for " + list(sLongestNuclei)[i])
            sheet.write(2, 3 + iSummaryErrors*len(sLongestNuclei) 
                            + i*iSummaryErrors, 
                        "Maximum Absolute Error for " 
                            + list(sLongestNuclei)[i])   
    
    if ShiftBar >= 1 and ShieldingBar >= 1: 
        iRightColumn = (len(sLongestNuclei)*iSummaryErrors 
                        + len(sLongestNucleiShift)*iSummaryErrors 
                        + len(sLongestNucleiShielding)*iSummaryErrors)
        if len(sLongestNucleiShielding) > 0:
            sheet.merge_range(1, len(sLongestNuclei)*iSummaryErrors 
                                  + len(sLongestNucleiShift)*iSummaryErrors + 1, 
                              1, len(sLongestNuclei)*iSummaryErrors 
                                  + len(sLongestNucleiShift)*iSummaryErrors + 1 
                                  + len(sLongestNucleiShielding)*iSummaryErrors - 1, 
                              "Errors of Scaled Shifts Shielding")  
        else:
            sheet.merge_range(1, len(sLongestNuclei)*iSummaryErrors 
                                  + len(sLongestNucleiShift)*iSummaryErrors + 1, 
                              1, len(sLongestNuclei)*iSummaryErrors 
                                  + len(sLongestNucleiShift)*iSummaryErrors + 1
                                  + len(sLongestNucleiShielding)*iSummaryErrors, 
                              "Errors of Scaled Shifts Shielding")
        for i in xrange(0,len(sLongestNucleiShielding)):
            sheet.write(2, 1 + len(sLongestNuclei)*iSummaryErrors 
                           + len(sLongestNucleiShift)*iSummaryErrors
                           + i*iSummaryErrors, 
                        "Mean Absolute Error for " 
                           + list(sLongestNuclei)[i])
            sheet.write(2, 2 + len(sLongestNuclei)*iSummaryErrors 
                           + len(sLongestNucleiShift)*iSummaryErrors 
                           + i*iSummaryErrors, 
                        "Root Mean Square Error for " 
                           + list(sLongestNuclei)[i])
            sheet.write(2, 3 + len(sLongestNuclei)*iSummaryErrors 
                           + len(sLongestNucleiShift)*iSummaryErrors 
                           + i*iSummaryErrors, 
                        "Maximum Absolute Error for " 
                           + list(sLongestNuclei)[i])  
           
    sheet.merge_range(0, 1, 0, iRightColumn, "Errors")
    sheet.write(2, iRightColumn + 1, "Average CPU times (s)")
    sheet.write(2, iRightColumn + 2, "Standard Deviation of CPU TImes (s)")
  
    iMoleculeNumber = SHEET1.nrows
    
    iSummaryErrors = 3
    iResultErrors = 5
    
    workbook2 = xlrd.open_workbook("Results.xlsx")
    
    for method_no in range(1, SHEET2.nrows):
    
        tOriginalParameters, tParameters = exporters.ExtractInputs(SHEET2, method_no)
        
        nrow = 3 + method_no - 1

        try: 
            sheet2 = workbook2.sheet_by_name("Method" + str(method_no)) 
        except:
            sheet2 = None
            print ("\nWARNING! No error calculation is performed for the Method "
                    + tools.method(tOriginalParameters).replace("_","/")) 
        
        if sheet2 != None:

            sheet.write(nrow, 0, str(tools.method(tParameters)))  
           
            for i in xrange(0, len(tParameters[4])):
                
                iTopRow = 3
                iBottomRow = 3 + iMoleculeNumber - 1
                #iBorromRow = 3+(sheet2.nrows-4)/2
                for j in xrange(1, len(sLongestNuclei)*iSummaryErrors+1): 
                    
                    inucleus = int(math.floor((j-1)/iSummaryErrors))
                    
                    if j % iSummaryErrors == 1:                   
                        if list(tParameters[4])[i] in sLongestNuclei[inucleus]:
                            sheet.write(nrow, j, 
                                        basicoperations.average_col(sheet2, iTopRow, 
                                                                    iBottomRow, 
                                                                    2+i*iResultErrors))
                    if j % iSummaryErrors == 2:
                        if list(tParameters[4])[i] in sLongestNuclei[inucleus]:
                            sheet.write(nrow, j, 
                                        basicoperations.average_col(sheet2, iTopRow, 
                                                                    iBottomRow, 
                                                                    4+i*iResultErrors))
                    if j % iSummaryErrors == 0:
                        if list(tParameters[4])[i] in sLongestNuclei[inucleus]:   
                            sheet.write(nrow, j, 
                                        basicoperations.average_col(sheet2, iTopRow, 
                                                                    iBottomRow, 
                                                                    6+i*iResultErrors))         
                
            if lParameters[method_no-1][0] == True: #if tParameters[7] == "YES":
                #if ShiftBar >= 1: 
                for i in xrange(0, len(tParameters[4])):
                    
                    iTopRow = 4 + iMoleculeNumber
                    #iBottomRow = sheet2.nrows/2+2
                    iBottomRow = 3 + 2*iMoleculeNumber 
                    #iBorromRow = sheet2.nrows
                    for j in xrange(len(sLongestNuclei)*iSummaryErrors + 1,
                                    len(sLongestNuclei)*iSummaryErrors 
                                        + len(sLongestNucleiShift)*iSummaryErrors + 1
                                    ):  
                        
                        inucleus = int(math.floor((j
                                            -len(sLongestNuclei)*iSummaryErrors-1)
                                        /iSummaryErrors))
                        
                        if j % iSummaryErrors == 1:
                            if list(tParameters[4])[i] in sLongestNucleiShift[inucleus]:
                                sheet.write(nrow, j, 
                                            basicoperations.average_col(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    2+i*iResultErrors))
                        if j % iSummaryErrors == 2:
                            if list(tParameters[4])[i] in sLongestNucleiShift[inucleus]:
                                sheet.write(nrow, j, 
                                            basicoperations.average_col(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    4+i*iResultErrors))
                        if j % iSummaryErrors == 0:
                            if list(tParameters[4])[i] in sLongestNucleiShift[inucleus]:   
                                sheet.write(nrow, j, 
                                            basicoperations.average_col(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    6+i*iResultErrors))

            if (lParameters[method_no-1][0] == False and 
                    lParameters[method_no-1][1] == True): 
                
                for i in xrange(0, len(tParameters[4])):

                    #if tParameters[7] != "YES" and tParameters[8] == "YES":            
                    iTopRow = 4 + iMoleculeNumber
                    #iBottomRow = sheet2.nrows/2+2
                    iBottomRow = 3 + 2*iMoleculeNumber 
                    #iBorromRow = sheet2.nrows
                    
                    for j in xrange(len(sLongestNuclei)*iSummaryErrors + 1, 
                                    len(sLongestNuclei)*iSummaryErrors 
                                      + len(sLongestNucleiShielding)*iSummaryErrors + 1
                                    ):  

                        inucleus = int(math.floor((j
                                            -len(sLongestNuclei)*iSummaryErrors-1)
                                        /iSummaryErrors))
                        
                        if j % iSummaryErrors == 1:
                            if list(tParameters[4])[i] in sLongestNucleiShielding[inucleus]:
                                sheet.write(nrow, j, 
                                            basicoperations.average_col(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    2+i*iResultErrors))
                        if j % iSummaryErrors == 2:
                            if list(tParameters[4])[i] in sLongestNucleiShielding[inucleus]:
                                sheet.write(nrow, j, 
                                            basicoperations.average_col(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    4+i*iResultErrors))
                        if j % iSummaryErrors == 0:
                            if list(tParameters[4])[i] in sLongestNucleiShielding[inucleus]:   
                                sheet.write(nrow, j, 
                                            basicoperations.average_col(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    6+i*iResultErrors))
                                
            if (lParameters[method_no-1][0] == True and 
                    lParameters[method_no-1][1] == True):  
                
                for i in xrange(0, len(tParameters[4])):
                    
                    #if tParameters[7] == "YES" and tParameters[8] == "YES":
                    iTopRow = 5 + 2*iMoleculeNumber
                    iBottomRow = 4 + 3*iMoleculeNumber 
                    
                    for j in xrange(len(sLongestNuclei)*iSummaryErrors 
                                      + len(sLongestNucleiShift)*iSummaryErrors + 1, 
                                    len(sLongestNuclei)*iSummaryErrors 
                                      + len(sLongestNucleiShift)*iSummaryErrors 
                                      + len(sLongestNucleiShielding)*iSummaryErrors + 1
                                    ):  
                        
                        inucleus = int(math.floor((j
                                            -len(sLongestNuclei)*iSummaryErrors
                                            -len(sLongestNucleiShift)*iSummaryErrors-1)
                                        /iSummaryErrors))
                        
                        if j % iSummaryErrors == 1:
                            if list(tParameters[4])[i] in sLongestNucleiShielding[inucleus]:
                                sheet.write(nrow, j, 
                                            basicoperations.average_col(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    2+i*iResultErrors))
                        if j % iSummaryErrors == 2:
                            if list(tParameters[4])[i] in sLongestNucleiShielding[inucleus]:
                                sheet.write(nrow, j, 
                                            basicoperations.average_col(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    4+i*iResultErrors))
                        if j % iSummaryErrors == 0:
                            if list(tParameters[4])[i] in sLongestNucleiShielding[inucleus]:   
                                sheet.write(nrow, j, 
                                            basicoperations.average_col(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    6+i*iResultErrors))
                
            iTopRow = 3
            iBottomRow = 3 + iMoleculeNumber - 1
            sheet.write(nrow, iRightColumn + 1, 
                                            basicoperations.average_col(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    7+i*iResultErrors))
            sheet.write(nrow, iRightColumn + 2, 
                                            basicoperations.std(sheet2, 
                                                                    iTopRow, 
                                                                    iBottomRow, 
                                                                    7+i*iResultErrors))

    workbook.close()

    print "\nSummaries.xlsx File has been Created!"

def error_file(SHEET1, SHEET2):
    
    lParameters = exporters.scalingparameters(SHEET2)
    
    workbook = xlsxwriter.Workbook('Results.xlsx')
    
    for method_no in range(1, SHEET2.nrows):
    
        tOriginalParameters, tParameters = exporters.ExtractInputs(SHEET2, method_no)

        sheet = workbook.add_worksheet("Method" + str(method_no))
        
        #merge_range(first_row, first_col, last_row, last_col, data[, cell_format])
        sheet.merge_range(0, 0, 0, 2+5*len(tParameters[4]), 
                          str(tools.method(tParameters).replace("/","_"))) # line 1
        
        sheet.write(1, 0, "InChI") # line 2, column 1
        sheet.write(1, 1, "InChIKey") # line 2, column 2
    
        sheet.merge_range(2, 0, 2, 1+5*len(tParameters[4]), 
                          "Calculated Shifts") # line 3         
        
        for i in range(0, len(tParameters[4])):
            sheet.write(1, 2+i*5, 
                        "Mean Absolute Error for " + list(tParameters[4])[i])
            sheet.write(1, 3+i*5, 
                        "Mean Square Error for " + list(tParameters[4])[i]) #line 2, column 4
            sheet.write(1, 4+i*5, 
                        "Root Mean Square Error for " + list(tParameters[4])[i]) #line 2, column 5
            sheet.write(1, 5+i*5, 
                        "Mean Signed Error for " + list(tParameters[4])[i]) #line 2, column 6
            sheet.write(1, 6+i*5, 
                        "Maximum Absolute Error for " + list(tParameters[4])[i]) #line 2, column 7        

        sheet.write(1, 6+(len(tParameters[4])-1)*5+1, "CPU time (s)")
        
        if tools.detectformat(tParameters[6]) == "InChI":
             refInChI = tParameters[6]    
             refInChIKey = inchi2input.inchi2inchikey(refInChI)
             sRefBaseName = refInChIKey
        else:
            sRefBaseName = tParameters[6].split(".xyz")[0] 
            
        if lParameters[method_no-1][0] == True: 
            sheet.merge_range(SHEET1.nrows + 3, 0, 
                              SHEET1.nrows + 3, 1+5*len(tParameters[4]), 
                              "Scaled Shifts") #line 4 + # of Molecules    
        if lParameters[method_no-1][0] == True and lParameters[method_no-1][1] == True: 
            sheet.merge_range(2 * SHEET1.nrows + 4, 0, 
                              2 * SHEET1.nrows + 4, 1+5*len(tParameters[4]), 
                              "Predicted Shifts by Linear Regression")  
        if lParameters[method_no-1][0] == False and lParameters[method_no-1][1] == True:
            sheet.merge_range(SHEET1.nrows + 3, 0, 
                              SHEET1.nrows + 3, 1+5*len(tParameters[4]), 
                              "Predicted Shifts by Linear Regression") #line 4 + # of Molecules           

        for row in range(0, SHEET1.nrows): #number of molecules = sheet.nrows - 1
            
            cell_input = SHEET1.cell_value(row, 0).encode('utf-8')
            inputtype = tools.detectformat(cell_input)
            if inputtype == "InChI":
                InChI, InChIKey = inchi2input.InchIs(SHEET1, row, 0)   
                sBaseName = InChIKey 
            elif inputtype == ".xyz":
                sBaseName = cell_input.split(".xyz")[0] 
                
            sFileName = (sBaseName + "_"
                            + str(tools.method(tParameters).replace("/","_")) 
                            + ".output")
        
            if tools.CheckExistLocal(sFileName) == False:
    
                print ("\nErrors could not be calculated for the molecule: "
                      + sBaseName + "\t" + tools.method(tOriginalParameters))
           
            elif tools.CheckText(sFileName, "Total times") == False:
                
                print ("\nErrors could not be calculated for the molecule: "
                      + sBaseName + "\t" + tools.method(tOriginalParameters))
            
            elif tools.CheckExistLocal(sBaseName + "_3D.mol") == False:
                
                print ("\nErrors could not be calculated for the molecule: "
                      + sBaseName + "\t" + tools.method(tOriginalParameters))
            else:
            
                if inputtype == "InChI": sheet.write(row + 3, 0, InChI)
                sheet.write(row + 3, 1, sBaseName)    
            
                if lParameters[method_no-1][0] == True:
                    if inputtype == "InChI": 
                        sheet.write(row + SHEET1.nrows + 4, 0, InChI)
                    sheet.write(row + SHEET1.nrows + 4, 1, sBaseName)  
                if lParameters[method_no-1][0] == True and lParameters[method_no-1][1] == True: 
                    if inputtype == "InChI": 
                        sheet.write(row + 2*SHEET1.nrows + 5, 0, InChI)
                    sheet.write(row + 2*SHEET1.nrows + 5, 1, sBaseName)
                if lParameters[method_no-1][0] == False and lParameters[method_no-1][1] == True: 
                    if inputtype == "InChI": 
                        sheet.write(row + SHEET1.nrows + 4, 0, InChI)
                    sheet.write(row + SHEET1.nrows + 4, 1, sBaseName) 
                  
                  
                for i in range(0, len(tParameters[4])):
                    
                    sBlock = "Calculated Shifts Relative to the Reference of: " \
                             + sRefBaseName + "\n" \
                             + tools.method(tOriginalParameters)
                    sheet.write(row + 3, 2+i*5, 
                                str(basicoperations.mae(sBaseName, sBlock, 
                                                        list(tParameters[4])[i])))
                    sheet.write(row + 3, 3+i*5, 
                                str(basicoperations.mse(sBaseName, sBlock, 
                                                        list(tParameters[4])[i])))            
                    sheet.write(row + 3, 4+i*5, 
                                str(basicoperations.rmse(sBaseName, sBlock, 
                                                         list(tParameters[4])[i])))
                    sheet.write(row + 3, 5+i*5, 
                                str(basicoperations.msie(sBaseName, sBlock, 
                                                         list(tParameters[4])[i])))           
                    sheet.write(row + 3, 6+i*5, 
                                str(basicoperations.maxae(sBaseName, sBlock, 
                                                          list(tParameters[4])[i])))               
                    
                    if lParameters[method_no-1][0] == True: 
                        
                        sBlock = "Scaled Shifts Relative to the Reference of: " \
                                 + sRefBaseName + "\n" \
                                 + tools.method(tOriginalParameters)
                        sheet.write(row + SHEET1.nrows + 4, 2+i*5, 
                                    str(basicoperations.mae(sBaseName, sBlock, 
                                                            list(tParameters[4])[i])))
                        sheet.write(row + SHEET1.nrows + 4, 3+i*5, 
                                    str(basicoperations.mse(sBaseName, sBlock, 
                                                            list(tParameters[4])[i])))          
                        sheet.write(row + SHEET1.nrows + 4, 4+i*5, 
                                    str(basicoperations.rmse(sBaseName, sBlock, 
                                                             list(tParameters[4])[i])))
                        sheet.write(row + SHEET1.nrows + 4, 5+i*5, 
                                    str(basicoperations.msie(sBaseName, sBlock, 
                                                             list(tParameters[4])[i])))           
                        sheet.write(row + SHEET1.nrows + 4, 6+i*5, 
                                    str(basicoperations.maxae(sBaseName, sBlock, 
                                                             list(tParameters[4])[i])))
                        
                    if lParameters[method_no-1][0] == True and lParameters[method_no-1][1] == True: 
                        
                        sBlock = "Predicted Shifts by Linear Regression"+ "\n" \
                                 + tools.method(tOriginalParameters)
                        sheet.write(row + 2*SHEET1.nrows + 5, 2+i*5, 
                                    str(basicoperations.mae(sBaseName, sBlock, 
                                                            list(tParameters[4])[i]))) 
                        sheet.write(row + 2*SHEET1.nrows + 5, 3+i*5, 
                                    str(basicoperations.mse(sBaseName, sBlock, 
                                                            list(tParameters[4])[i])))
                        sheet.write(row + 2*SHEET1.nrows + 5, 4+i*5, 
                                    str(basicoperations.rmse(sBaseName, sBlock, 
                                                             list(tParameters[4])[i])))
                        sheet.write(row + 2*SHEET1.nrows + 5, 5+i*5, 
                                    str(basicoperations.msie(sBaseName, sBlock, 
                                                             list(tParameters[4])[i])))           
                        sheet.write(row + 2*SHEET1.nrows + 5, 6+i*5, 
                                    str(basicoperations.maxae(sBaseName, sBlock, 
                                                             list(tParameters[4])[i])))
                        
                    if lParameters[method_no-1][0] == False and lParameters[method_no-1][1] == True: 
                        
                        sBlock = "Predicted Shifts by Linear Regression"+ "\n" \
                                 + tools.method(tOriginalParameters)
                        sheet.write(row + SHEET1.nrows + 4, 2+i*5, 
                                    str(basicoperations.mae(sBaseName, sBlock, 
                                                            list(tParameters[4])[i])))
                        sheet.write(row + SHEET1.nrows + 4, 3+i*5, 
                                    str(basicoperations.mse(sBaseName, sBlock, 
                                                            list(tParameters[4])[i])))          
                        sheet.write(row + SHEET1.nrows + 4, 4+i*5, 
                                    str(basicoperations.rmse(sBaseName, sBlock, 
                                                             list(tParameters[4])[i])))
                        sheet.write(row + SHEET1.nrows + 4, 5+i*5, 
                                    str(basicoperations.msie(sBaseName, sBlock, 
                                                             list(tParameters[4])[i])))           
                        sheet.write(row + SHEET1.nrows + 4, 6+i*5, 
                                    str(basicoperations.maxae(sBaseName, sBlock, 
                                                             list(tParameters[4])[i])))
           
                sCPUtime = exporters.CPUtime(sFileName.split(".output")[0])
                sCPUtime = (sCPUtime.split("s")[0] 
                            if "s" in sCPUtime
                            else sCPUtime)
                sheet.write(row + 3, 6+(len(tParameters[4])-1)*5+1, sCPUtime)
           
    workbook.close() 
    
    print "\nResults.xlsx File has been Created!" 
 
def AddExpState(sBaseName):    
    '''
    It writes the experimental specs to the 3D mol file of 
    a given molecule
    '''
    
    with open(sBaseName + "_3D.mol", 'a+') as file:
       file.write("\nExperimental Shifts"
                 + "\nFT spectra: \tSpeed: \tSpectra resoultion: " 
                 + "\tFrequency: \tSolvent: \tTemperature: \tReference: \n")
       [file.write(str(exporters.NucleiNoandListXYZ(sBaseName)[i][0]) + "\t" + 
                   str(exporters.NucleiNoandListXYZ(sBaseName)[i][1]) + "\t\n") 
                   for i in xrange(0, 
                                  len(exporters.NucleiNoandListXYZ(sBaseName)))
                   ]
 
def CreateNWCHEMFile(sBaseName, sLocal, sSource, tOriginalParameters, tParameters):
    
    ''' It prepares the .nw file, an input file submitted to NWChem'''
    
    filename = (sBaseName + "_" 
                + str(tools.method(tParameters).replace("/","_")) 
                + ".nw")

    # solvent
    if tOriginalParameters[5] == 'GAS':         
        shutil.copyfile(sLocal + "/isicle-packages/temp_gas.txt", sLocal + '//' + filename)        
    else:        
        shutil.copyfile(sLocal + "/isicle-packages/temp_solvent.txt", sLocal + '//' + filename)        
        tools.ReplaceText(filename, "SSS", tOriginalParameters[5])
    
    tools.ReplaceText(filename, "MMM",
                      sBaseName + "_" 
                      + str(tools.method(tParameters).replace("/","_")))
    
    # xyz file input & driver
    tools.ReplaceText(filename, "sSource", sSource)   
    
    iNwchemCharLimit = 60
    sNwchemGeoFileName = sBaseName + tParameters[2].upper() + tParameters[3].upper()
    
    if len(sNwchemGeoFileName) > iNwchemCharLimit:        
        tools.ReplaceText(filename, "XXX_GGG1_GGG2", sNwchemGeoFileName[:60])          
    else:
        # xyz file input & driver
        tools.ReplaceText(filename, "XXX", sBaseName)
    
        # xyz file input & driver # geometry functional
        tools.ReplaceText(filename, "GGG1", tParameters[2].upper())
    
        # xyz file input & driver # geometry basis
        tools.ReplaceText(filename, "GGG2", tParameters[3].upper()) 
    
    # xyz file input & driver
    tools.ReplaceText(filename, "XXX", sBaseName)
        
    # geometry basis
    tools.ReplaceText(filename, "YYY1", tOriginalParameters[3]) 
    
    # geometry functional
    tools.ReplaceText(filename, "ZZZ1", tOriginalParameters[2]) 
    
    # nmr basis
    tools.ReplaceText(filename, "YYY2", tOriginalParameters[1]) 
    
    # nmr functional
    tools.ReplaceText(filename, "ZZZ2", tOriginalParameters[0]) 
    
    tools.ReplaceText(filename, "PPP", 
                      exporters.shielding_index(sBaseName, 
                                                tOriginalParameters[4]))

def print_methods(sheet2):
    
    for i in range(1, sheet2.nrows):
        
        tOriginalParameters, tParameters  = exporters.ExtractInputs(sheet2, i)
    
        print ("\n" + "-"*79 + "\n" + "Method No: " + str(i)  
           + " " + tools.method(tOriginalParameters))     
        