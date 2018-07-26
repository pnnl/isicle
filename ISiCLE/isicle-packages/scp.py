# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 23:35:58 2016

@author: yesi172
"""

import getpass
import sys
import os
import time
import shutil
import hashlib

from fabric.api import env
from fabric.operations import run, put, get, hide, settings

import tools
import inchi2input
import basicoperations
import exporters

def get_userinfo(sHOST):

    '''This function gets user id and password'''

    env.host_string = sHOST
    env.user        = raw_input("Network ID: ")
    env.password    = getpass.getpass('Password: ')

def clean_userinfo():

    '''This function deletes user id and password'''

    env.host_string = ""
    env.user        = ""
    env.password    = ""

def upload_files(sLocalSource, sServerDir):

    '''
    This function copies a given folder or file in local computer
    to server.  sLocalSource: Directory or a File name in local
    computer.  sServerDir: Directory where the given folder or
    file is copied in server
    '''
    try:        
        with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                      warn_only=True):
            #run('mkdir -p ' + sDest)  #make sure the directory is there!
            #put(sLocalSource, sServerDir, mode=0664).succeeded
            put(sLocalSource, sServerDir)
    except:
        print sLocalSource + " could not be uploaded!"
        
def download_files(sServerSource, sLocalDest):

    '''
    This function copies a given folder or file in server to local
    computer.  sServerSource: Directory or a File name in server
    sLocalDest: Directory where the given folder or file is copied
    in local computer.
    '''

    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                      warn_only=True):

        #local('mkdir -p ' +  sLocalDest + '/scpTest/test')

        get(sServerSource, sLocalDest)

def _run(sHOST, sServerDir):

    '''
    This function submits job as stated in the file named setrun.bash
    on Cascade.
    '''

    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                  warn_only=True):
        env.warn_only = True
        run("tr -d '\r' < " + sServerDir + "/setrun.bash > "
            + sServerDir + "/finalsetrun.bash")
        #run("cat "+ sDest + "/setrun.bash  | tr -d '\\r'  > " + sDest +"/test" + "/finalsetrun.bash")
        #tr = trancate #-d = delete new line from setrun.bash
        #\r = new line #> means write them into newfile  #cat =
        
    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                  warn_only=True):        
        run("chmod +x " + sServerDir + "finalsetrun.bash") #changemod #+x means executable file
        
    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                  warn_only=True):        
        jobInfo = run("msub " + sServerDir + "finalsetrun.bash")

    if "Submitted batch job " in jobInfo:
            jobID = jobInfo.split("Submitted batch job ")[-1]
            
    return jobID

def __run(sHOST, sServerDir, sFileName):

    '''This function submits a job on Cascade'''

    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                      warn_only=True):
        env.warn_only = True

        jobID = run("submit_nwchem /" + sServerDir + "/" + sFileName
                    + ".nw -accnt emslc49661 -nn 1 -nc 16 -time 1:00:00:00")
    
    return jobID

def check_job_status(jobID):
    
    jobInfo = jobstatus(jobID)
    if type(jobInfo) == tuple and jobInfo[2] == "Deferred": 
        # Blocked job
        sys.exit("NOTE: Cascade resource manager is currently down")
        
    #while jobstatus(jobID)[2] == "Idle":
    #    time.sleep(10)
    
    jobInfo = jobstatus(jobID)  
    bReportJobOnce = False
    while type(jobInfo) == tuple and jobInfo[2] == "PD":
        if bReportJobOnce == False:
            print jobInfo[0] 
            bReportJobOnce = True
        time.sleep(60)
        jobInfo = jobstatus(jobID)  
        
    #if jobstatus(jobID)[2] == "Running":
    #    print ("\nSimulation Has Been Started with Job ID of "
    #          + jobID + " on Cascade!\n")    
        
    jobInfo = jobstatus(jobID)        
    if type(jobstatus(jobID)) == tuple and jobInfo[2] == "R":
        print ("\nSimulation Has Been Started with Job ID of "
              + jobID + " on Cascade!\n")  
   
    #jobInfo = jobstatus(jobID)        
    #if (type(jobInfo) == tuple and 
    #   (jobInfo[2] != "R" and jobInfo[2] != "PD")):
    #    print "NOTE: No connection is available to service this operation"
        
def jobstatus(jobID):

    '''
    This functions checks the status of a submitted job on Cascade
    by its job ID
    Command:  showq | grep jobID
    sJoInfo: output of the command above
    1st element: Job ID
    2nd element: Username
    3rd element: Number of nodes used -> 1, 2, 16
    4th element: JobStatus -> Idle, Running, PD, R, ...
    5th element: Remanining time -> 00:10:00, ....
    6th element:Starting time -> mm-ddDhh:mm:ss
    '''

    try:
        with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                      warn_only=True):
            sJobInfo = run("showq | grep " + jobID)
            #if 'Option w requires an argument' in sJobInfo:
            #    lJobInfo = sJobInfo.split('Option w requires an argument\r\n')[-1] #List

            lJobInfo = sJobInfo.split()
            if len(lJobInfo) >= 6: #if jobID is active
                sJobID = lJobInfo[0]
                sUserName = lJobInfo[1]
                sNodeNumber = lJobInfo[2]
                sJobState = lJobInfo[3]
                sJobRemaining = lJobInfo[4]
                return sJobInfo, sUserName, sJobState, sNodeNumber, sJobRemaining
            else:
                print jobID, "is not active"
                return "" #if jobID is not active
    except:
        return "" 
        
def prepare_runfile(sACCOUNT, iNODES, iPROCESSORS, sWalltime, sDest, sLocalDir):

    '''
    This function takes the number of nodes, number of processors,
    wall time given by user and writes them to the file named
    setrun.bash which will be used to submit the jobs.
    '''

    shutil.copyfile(sLocalDir + "/isicle-packages/temp_setrun.bash", sLocalDir + "/isicle-packages/setrun.bash")

    tools.ReplaceText(sLocalDir + "/isicle-packages/setrun.bash", 
                      "#MSUB -A sAccount",
                      "#MSUB -A " + sACCOUNT)
    tools.ReplaceText(sLocalDir + "/isicle-packages/setrun.bash", 
                      "#MSUB -l nodes=__:ppn=__,walltime=__",
                      "#MSUB -l nodes=" + iNODES + ":ppn="
                          + iPROCESSORS + ",walltime=" + sWalltime)
    tools.ReplaceText(sLocalDir + "/isicle-packages/setrun.bash", 
                      "#MSUB -o __.%j.out",
                      "#MSUB -o " + sDest + ".%j.out")
    tools.ReplaceText(sLocalDir + "/isicle-packages/setrun.bash", 
                      "#MSUB -e __.%j.err",
                      "#MSUB -e " + sDest + ".%j.err")
    tools.ReplaceText(sLocalDir + "/isicle-packages/setrun.bash", 
                      "-n __", "-n " + iPROCESSORS)
    tools.ReplaceText(sLocalDir + "/isicle-packages/setrun.bash", 
                      "cd __", "cd " + sDest)
        
    shutil.copyfile(sLocalDir + "/isicle-packages/setrun.bash", sLocalDir + "/setrun.bash")
    
def get_walltime():

    '''
    This function takes a new time wall time from user if the current
    job stops running

    '''
    print ("\nError: Wall time has been reached or job has become "
            + "unsuccessful before simulation is completed")    
    
    sWallTime = ""
    bCont = False
    
    while sWallTime is "":

        UserChoice = raw_input("\nWould you like to continue simulation (Y/N): ")
        if UserChoice.upper() == "Y":
            bCont = True
            sWallTime = raw_input("\nEnter the wall time (dd:hh:mm:ss): ")
            if tools.detect_walltime_format(sWallTime) == False:
                sWallTime = ""
                print "\nWall time was not given in reqired format!"
        elif UserChoice.upper() == "N":
            print "\nSimulation is terminated by user!"
            bCont = False
            return bCont, ""
        else:
            print "\nInvalid Character"
            continue

    return bCont, sWallTime

def check_text_server(sFileName, sText):

    '''
    Function to check if a given text is found in the name of the
    file, sFileName.
    '''

    try:
        with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                                  warn_only=True):
            return run('grep ' + sText + ' ' + sFileName)
    except:
        return ""

def globserver(sFileNameServer):

    '''This fucntion lists the files beginning with sFileName on server'''

    try:
        with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                      warn_only=True):
            return run('ls ' + sFileNameServer) #FileList
    except:
        return ""

def replace_invalidcharacters(sFileName):

    '''
    This function replaces the charachter of asterisk *, which is a
    reserved character on Windows operating system
    '''

    if "*" in sFileName:
        FileList = globserver(sFileName)
        FileList = FileList.split()
        if FileList != "":
            for sFileName in FileList:
                sNewFileName = sFileName.replace("*","STAR")
                with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                                  warn_only=True):
                    run("mv " + sFileName + " " + sNewFileName)
                    #run("rename 's/\*/+/g' " + sFileName + " -vn") #FileList

def check_filessame(sServerFileName, sLocalFileName):

    '''
    This function checks if two files in server and local computer are
    same by using hashing method
    '''

    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                      warn_only=True):

        if (run("md5sum " + sServerFileName).split()[0]
            == hashlib.md5(open(sLocalFileName, 'rb').read()).hexdigest()):
            return True
        else:
            return False
            

def check_filesdownloaded(sServerSource, sLocalDest, sheet1, sheet2):
    '''
    This function checks that the output files are downloaded from server
    to local computer.  It ckeck the files with the extension of .output,
    _3D mol and .xls.
    '''

    # Output files
    for i in range(1,sheet2.nrows):
        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i)
        if tools.detectformat(tParameters[6]) == "InChI":
            refInChI = tParameters[6]
            refInChIKey = inchi2input.inchi2inchikey(refInChI)
            sRefBaseName = refInChIKey
        else:
            sRefBaseName = tParameters[6].split(".xyz")[0]
        CompletedCount = 0
        # Reference molecule
        sFileNameLocal = sLocalDest + "/" + sRefBaseName + "_"\
                            + str(tools.method(tParameters).replace("/","_"))\
                            + ".output"
        sFileNameServer = sServerSource + '/' + sRefBaseName + '_'\
                            + tools.method(tParameters).replace("/","_")\
                            + ".output"
        if redownload(sFileNameLocal, sFileNameServer,
                      sLocalDest, sServerSource) == True:
            CompletedCount = CompletedCount + 1
        # Molecules
        if CompletedCount <= sheet1.nrows:
            #while CompletedCount <= sheet1.nrows-1: +1 comes from reference
            for j in range(0, sheet1.nrows):
                cell_input = sheet1.cell_value(j, 0).encode('utf-8')
                if tools.detectformat(cell_input) == "InChI":
                    InChI = cell_input
                    InChIKey = inchi2input.inchi2inchikey(InChI)
                    sBaseName = InChIKey
                else:
                    sBaseName = cell_input.split(".xyz")[0]
                sFileNameLocal = (sLocalDest + "/" + sBaseName + "_"
                                 + str(tools.method(tParameters).replace("/","_"))
                                 + ".output")
                sFileNameServer = (sServerSource + "/" + sBaseName + '_'
                                  + tools.method(tParameters).replace("/","_")
                                  + ".output")
                if redownload(sFileNameLocal, sFileNameServer,
                              sLocalDest, sServerSource) == True:
                    CompletedCount = CompletedCount + 1

    # 3D Mol files
    for i in range(1, sheet2.nrows):
        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i)
        if tools.detectformat(tParameters[6]) == "InChI":
            refInChI = tParameters[6]
            refInChIKey = inchi2input.inchi2inchikey(refInChI)
            sRefBaseName = refInChIKey
        else:
            sRefBaseName = tParameters[6].split(".xyz")[0]
        # Reference molecule
        CompletedCount = 0
        sFileNameLocal = sLocalDest + "/" + sRefBaseName + "_3D.mol"
        sFileNameServer = sServerSource + "/" + sRefBaseName + "_3D.mol"
        if redownload(sFileNameLocal, sFileNameServer,
                      sLocalDest, sServerSource) == True:
            CompletedCount = CompletedCount + 1
        # Molecules
        if CompletedCount <= sheet1.nrows:
            #while CompletedCount <= sheet1.nrows-1: +1 comes from reference
            for j in range(0, sheet1.nrows):
                cell_input = sheet1.cell_value(j, 0).encode('utf-8')
                if tools.detectformat(cell_input) == "InChI":
                    InChI = cell_input
                    InChIKey = inchi2input.inchi2inchikey(InChI)
                    sBaseName = InChIKey
                else:
                    sBaseName = cell_input.split(".xyz")[0]
                sFileNameLocal = sLocalDest + "/" + sBaseName + "_3D.mol"
                sFileNameServer = sServerSource + '/' + sBaseName + "_3D.mol"
                if redownload(sFileNameLocal, sFileNameServer,
                              sLocalDest, sServerSource) == True:
                              CompletedCount = CompletedCount + 1

    # 3D Mol files
    for i in range(1, sheet2.nrows):
        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i)
        if tools.detectformat(tParameters[6]) == "InChI":
            refInChI = tParameters[6]
            refInChIKey = inchi2input.inchi2inchikey(refInChI)
            sRefBaseName = refInChIKey
        else:
            sRefBaseName = tParameters[6].split(".xyz")[0]
        # Reference molecule
        CompletedCount = 0
        sFileNameLocal = sLocalDest + "/" + sRefBaseName + "_3D.mol"
        sFileNameServer = sServerSource + "/" + sRefBaseName + "_3D.mol"
        if redownload(sFileNameLocal, sFileNameServer,
                      sLocalDest, sServerSource) == True:
            CompletedCount = CompletedCount + 1
        # Molecules
        if CompletedCount <= sheet1.nrows:
            #while CompletedCount <= sheet1.nrows-1: +1 comes from reference
            for j in range(0, sheet1.nrows):
                cell_input = sheet1.cell_value(j, 0).encode('utf-8')
                if tools.detectformat(cell_input) == "InChI":
                    InChI = cell_input
                    InChIKey = inchi2input.inchi2inchikey(InChI)
                    sBaseName = InChIKey
                else:
                    sBaseName = cell_input.split(".xyz")[0]
                sFileNameLocal = sLocalDest + "/" + sBaseName + "_3D.mol"
                sFileNameServer = sServerSource + '/' + sBaseName + "_3D.mol"
                if redownload(sFileNameLocal, sFileNameServer,
                              sLocalDest, sServerSource) == True:
                              CompletedCount = CompletedCount + 1

def report_generatedfiles_server(sServerDir, sLocalDir,
                                 sheet1, sheet2, jobID):

    '''
    This function reports the molecules whose input files are
    generated on Cascade.
    '''

    for i in range(1,sheet2.nrows):

        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i)

        if tools.detectformat(tParameters[6]) == "InChI":
            refInChI = tParameters[6]
            refInChIKey = inchi2input.inchi2inchikey(refInChI)
            sRefBaseName = refInChIKey
        else:
            sRefBaseName = tParameters[6].split(".xyz")[0]

        if globserver(sServerDir + sLocalDir.split("/")[-2] + "/"
                     + sRefBaseName + "_"
                     + str(tools.method(tParameters).replace("/","_"))
                     + ".output") != "":
            print ("\nReference Molecule Has Been Started!" + "\t"
                  + tParameters[6] + "\t"
                  + str(tools.method(tOriginalParameters)) + "\n")
        else:
            while (type(jobstatus(jobID)) == tuple and
                  basicoperations.walltime_str2int(jobstatus(jobID)[-1])
                      >= basicoperations.walltime_str2int("0:00:00:00")):
                if globserver(sServerDir + sLocalDir.split("/")[-2] + "/"
                             + sRefBaseName + "_"
                             + str(tools.method(tParameters).replace("/","_"))
                             + ".output") == "":
                    time.sleep(10)
                else:
                    print ("\nReference Molecule Has Been Started!" + "\t"
                          + tParameters[6] + "\t"
                          + str(tools.method(tOriginalParameters)) + "\n")
                    continue

        for irow in xrange(0, sheet1.nrows):
            # number of molecules = sheet.nrows - 1
            cell_input = sheet1.cell_value(irow, 0).encode('utf-8')
            if tools.detectformat(cell_input) == "InChI":
                InChI = cell_input
                InChIKey = inchi2input.inchi2inchikey(InChI)
                sBaseName = InChIKey
            else:
                sBaseName = cell_input.split(".xyz")[0]

            if globserver(sServerDir + sLocalDir.split("/")[-2] + "/"
                         + sBaseName + "_"
                         + str(tools.method(tParameters).replace("/","_"))
                         + ".output") != "":
                imolecule = irow + 1
                print ("\nMolecule " + str(imolecule)
                      + " Has Been Started!" + "\t" + cell_input + "\t"
                      + str(tools.method(tOriginalParameters)) + "\n")
            else:
                while (type(jobstatus(jobID)) == tuple and
                      basicoperations.walltime_str2int(jobstatus(jobID)[-1])
                          >= basicoperations.walltime_str2int("0:00:00:00")):
                    if globserver(sServerDir + sLocalDir.split("/")[-2] + "/"
                                 + sBaseName + "_"
                                 + str(tools.method(tParameters).replace("/","_"))
                                 + ".output*") == "":
                        time.sleep(10)
                        irow -= 1
                    else:
                        imolecule = irow + 1
                        print ("\nMolecule " + str(imolecule)
                              + " Has Been Started!" + "\t" + cell_input + "\t"
                              + str(tools.method(tOriginalParameters)) + "\n")

def remove_incompletefiles_server(sLocalDir, sServerDir, sheet1, sheet2):

    for i in range(1, sheet2.nrows):

        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i)

        if tools.detectformat(tParameters[6]) == "InChI":
            refInChI = tParameters[6]
            refInChIKey = inchi2input.inchi2inchikey(refInChI)
            sRefBaseName = refInChIKey
        else:
            sRefBaseName = tParameters[6].split(".xyz")[0]
        # Reference
        remove_incomp_server(sServerDir + sLocalDir.split("/")[-2] + "/"
                            + sRefBaseName + "_"
                            + str(tools.method(tParameters).replace("/","_"))
                            + ".output")

        for irow in range(0, sheet1.nrows):
            cell_input = sheet1.cell_value(irow, 0).encode('utf-8')
            if tools.detectformat(cell_input) == "InChI":
                InChI = cell_input
                InChIKey = inchi2input.inchi2inchikey(InChI)
                sBaseName = InChIKey
            else:
                sBaseName = cell_input.split(".xyz")[0]
            remove_incomp_server(sServerDir + sLocalDir.split("/")[-2] + "/"
                                + sBaseName + "_"
                                + str(tools.method(tParameters).replace("/","_"))
                                + ".output")

def remove_incomp_server(sFileNameServer):

    '''This function deletes the incompleted output files on server'''

    FileList = globserver(sFileNameServer)
    FileList = FileList.split()

    if FileList != "":
        for FileName in FileList:
            with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                          warn_only=True):
                output = run("grep Total " + FileName + " | grep times")
                if output == "":                    
                    with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                          warn_only=True):
                        run('rm ' + FileName)

def check_outputcompleted_server(sFileNameServer):

    '''
    This function checks if the output file for a given molecule is
    completed on server.
    '''

    bCIsValid = True #means done!

    FileList = globserver(sFileNameServer)
    
    if "No such file or directory" in FileList:
        bCIsValid = False 
      
    FileList = FileList.split() 

    if FileList != "":
        for FileName in FileList:
            try:
                with settings(hide('warnings', 'running', 'stdout', 'stderr'),
                          warn_only=True):
                    output = run("grep Total " + FileName + " | grep times")
                    if output == "":
                        bCIsValid = False
            except:
                bCIsValid = False

    return bCIsValid

def redownload(sFileNameLocal, sFileNameServer, sLocalDest, sServerSource):
    """
    This function checks the given files are completelety downloaded
    or not.  If not, it raises a warning message.
    """
    # First download it!
    if tools.CheckExistLocal(sFileNameLocal) == False:
        download_files(sFileNameServer, sLocalDest)
    chanceToDownload = 0
    while (chanceToDownload < 10 and 
          tools.CheckExistLocal(sFileNameLocal) == False):
        time.sleep(5)
        chanceToDownload = chanceToDownload + 1
    if (chanceToDownload == 10 and
       tools.CheckExistLocal(sFileNameLocal) == False):
        print "WARNING! " + sFileNameServer + " file could not be downloaded!"
    # Check that the files are same if downladed!
    # This step is for the cases when files are partially downladed
    chanceToDownload = 0
    if tools.CheckExistLocal(sFileNameLocal) == True:
        while (chanceToDownload < 10 and
              check_filessame(sFileNameServer, sFileNameLocal) == False):
            time.sleep(5)
            chanceToDownload = chanceToDownload + 1
        if check_filessame(sFileNameServer, sFileNameLocal) == True:
            return True
        else:
            print ("WARNING! " + sFileNameServer
                  + " file could not be downloaded!")
            return False

def resimulate(sheet1, sheet2, sServerDir, sLocalDir, sHOST, 
               iNodesNo, iProcessorNoPerNode, chanceToNewRun):

    bReportOnce = False
    bReportJobOnce = False
    for i in range(1,sheet2.nrows):

        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i)

        if tools.detectformat(tParameters[6]) == "InChI":
            refInChI = tParameters[6]
            refInChIKey = inchi2input.inchi2inchikey(refInChI)
            sRefBaseName = refInChIKey
        else:
            sRefBaseName = tParameters[6].split(".xyz")[0]
            
        if (check_outputcompleted_server(sServerDir + sLocalDir.split("/")[-2] 
                   + "/" + sRefBaseName + "_"
                   + str(tools.method(tParameters).replace("/","_"))
                   + ".output") == False and
               chanceToNewRun < 3):  # reference            
            bCont, sWalltime = get_walltime()            
            if bCont == True: 
                remove_incompletefiles_server(sLocalDir, sServerDir, sheet1, sheet2)
                prepare_runfile(str(iNodesNo), str(iProcessorNoPerNode), sWalltime,
                                sServerDir + sLocalDir.split("/")[-2] + "/", sLocalDir)
                upload_files(sLocalDir, sServerDir)
                jobID = _run(sHOST, sServerDir + sLocalDir.split("/")[-2] + "/")
                print ("\nJob Has Been Submitted to Cascade "
                      + "with the Job ID of " + jobID)
                check_job_status(jobID)
                if bReportOnce == False:
                    report_generatedfiles_server(sServerDir, sLocalDir, 
                                                 sheet1, sheet2, jobID)
                    bReportOnce = True
                jobInfo = jobstatus(jobID)  
                while type(jobInfo) == tuple: 
                    # while jobstatus(jobID)[2] == "Running"
                    if bReportJobOnce == False:
                        print jobInfo[0]
                        bReportJobOnce = True
                    time.sleep(60)
                    jobInfo = jobstatus(jobID) 
                chanceToNewRun = chanceToNewRun + 1
            else:
                chanceToNewRun = 3
                
    bReportOnce = False
    bReportJobOnce = False
    for i in range(1,sheet2.nrows):

        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i)

        for irow in xrange(0,sheet1.nrows):

            cell_input = sheet1.cell_value(irow, 0).encode('utf-8')
            if tools.detectformat(cell_input) == "InChI":
                InChI = cell_input
                InChIKey = inchi2input.inchi2inchikey(InChI)
                sBaseName = InChIKey
            else:
                sBaseName = cell_input.split(".xyz")[0]

            if (check_outputcompleted_server(sServerDir 
                    + sLocalDir.split("/")[-2] + "/" + sBaseName + "_"
                    + str(tools.method(tParameters).replace("/","_"))
                    + ".output") == False and
                chanceToNewRun < 3):
                bCont, sWalltime = get_walltime()
                if bCont == True: 
                    remove_incompletefiles_server(sLocalDir, sServerDir, sheet1, sheet2)
                    prepare_runfile(str(iNodesNo), str(iProcessorNoPerNode),
                                    sWalltime,
                                    sServerDir + sLocalDir.split("/")[-2] + "/", 
                                    sLocalDir)
                    upload_files(sLocalDir, sServerDir)
                    jobID = _run(sHOST, sServerDir + sLocalDir.split("/")[-2] + "/")
                    print ("\nJob Has Been Submitted to Cascade "
                          + "with the Job ID of " + jobID)
                    check_job_status(jobID)
                    if bReportOnce == False:
                        report_generatedfiles_server(sServerDir, sLocalDir, 
                                                     sheet1, sheet2, jobID)
                        bReportOnce = True
                    jobInfo = jobstatus(jobID)  
                    while type(jobInfo) == tuple: 
                        # while jobstatus(jobID)[2] == "Running"
                        if bReportJobOnce == False:
                            print jobInfo[0]
                            bReportJobOnce = True
                        time.sleep(60)    
                        jobInfo = jobstatus(jobID) 
                    chanceToNewRun = chanceToNewRun + 1
                else:
                    chanceToNewRun = 3
                    
    return (sServerDir, sLocalDir, sHOST,
            iNodesNo, iProcessorNoPerNode, chanceToNewRun)

def simulate(sheet1, sheet2, sLocalDir, sACCOUNT, sHOST, sServerDir,
             iNodesNo, iProcessorNoPerNode, sWalltime):
 
    #sheet1: Molecule Set, #sheet2: Methods
    chanceToNewRun = 0
    inchi2input.CreateAllInputs(sheet1, sheet2, sLocalDir, sServerDir,
                                            iNodesNo, iProcessorNoPerNode, sWalltime)

    prepare_runfile(sACCOUNT, str(iNodesNo), str(iProcessorNoPerNode), sWalltime,
                    sServerDir + sLocalDir.split("/")[-2] + "/", sLocalDir)

    upload_files(sLocalDir, sServerDir)
    print "\nFiles Have Been Uploaded to Cascade!"
   
    #tools.RemoveInCompleteOutputFilesLocal(sheet1, sheet2, sLocalDir)
    #remove_incompletefiles_server(sLocalDir, sServerDir, sheet1, sheet2)

    jobID = _run(sHOST, sServerDir + sLocalDir.split("/")[-2] + "/")
    print "\nJob Has Been Submitted to Cascade with the Job ID of " + jobID
   
    # Monitor the job
    check_job_status(jobID)
    jobInfo = jobstatus(jobID)  
    bReportJobOnce = False
    while type(jobInfo) == tuple: 
        # while jobstatus(jobID)[2] == "Running"
        if bReportJobOnce == False:
            report_generatedfiles_server(sServerDir, sLocalDir,
                                 sheet1, sheet2, jobID)
            print jobInfo[0]
            bReportJobOnce = True
        time.sleep(60)
        jobInfo = jobstatus(jobID) 
    if type(jobstatus(jobID)) != tuple:
        print ("Simulation has been completed.")

    chanceToNewRun = 1
    while chanceToNewRun < 3:
        sServerDir, sLocalDir, sHOST, \
        iNodesNo, iProcessorNoPerNode, chanceToNewRun \
            = resimulate(sheet1, sheet2, sServerDir, sLocalDir, sHOST,
                         iNodesNo, iProcessorNoPerNode, chanceToNewRun)
        chanceToNewRun = chanceToNewRun + 1
   
    for i in range(1, sheet2.nrows):
        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i)
        print ("\nSimulation Has Been Completed for the Method  " + str(i)
              + ": " + tools.method(tOriginalParameters) + "\n")

    for i in range(1, sheet2.nrows):

        tOriginalParameters, tParameters = exporters.ExtractInputs(sheet2, i)

        if tools.detectformat(tParameters[6]) == "InChI":
            refInChI = tParameters[6]
            refInChIKey = inchi2input.inchi2inchikey(refInChI)
            sRefBaseName = refInChIKey
        else:
            sRefBaseName = tParameters[6].split(".xyz")[0]
        # reference
        replace_invalidcharacters(sServerDir + sLocalDir.split("/")[-2] + "/"
            + sRefBaseName + "_"
            + str(tools.method(tOriginalParameters).replace("/","_")))

        for irow in xrange(0, sheet1.nrows):
            cell_input = sheet1.cell_value(irow, 0).encode('utf-8')
            if tools.detectformat(cell_input) == "InChI":
                InChI = cell_input
                InChIKey = inchi2input.inchi2inchikey(InChI)
                sBaseName = InChIKey
            else:
                sBaseName = cell_input.split(".xyz")[0]
            replace_invalidcharacters(sServerDir + sLocalDir.split("/")[-2]
                + "/" + sBaseName + "_"
                + str(tools.method(tOriginalParameters).replace("/","_")))
  
    #download_files(sDest + sLocalDir.split("/")[-2] + "/", sLocalDir + "../")
    sFolderName = sLocalDir + "/Results_" \
        + tools.DiscardCharacter(time.strftime("%d/%m/%Y"),"/") \
        + "_" + tools.DiscardCharacter(time.strftime("%H:%M:%S"),":") + "/"

    if not os.path.exists(sFolderName):
        os.makedirs(sFolderName)
    download_files(sServerDir + sLocalDir.split("/")[-2] + "/*", sFolderName)
    
    check_filesdownloaded(sServerDir + sLocalDir.split("/")[-2],
                          sFolderName, sheet1, sheet2)
       
    print "\nFiles Have Been Downloaded from Cascade!\n"
    os.chdir(sFolderName.split("/")[-2])
