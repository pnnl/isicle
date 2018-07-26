# -*- coding: utf-8 -*-
"""
Created on Wed May 25 16:09:23 2016

@author: yesi172
"""
import argparse
import time
import sys
import os
import site 

if __file__:
    sys.path.insert(0,
                    os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), 
                                                 "isicle-packages"))  
else:
    site.addsitedir(os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                                                 "isicle-packages"))

import exporters
import tools
import inchi2input
import loggers
import scp
import basicoperations

Start = time.time()
parser = argparse.ArgumentParser(description='Short sample app')

#parser.add_argument('-dir', action = "store", dest = "sLocalDir", 
#                    help = "File directory in local computer,"
#                           + " Ex: -dir C:/Users/UserName/Desktop/test/")
parser.add_argument('-accnt', action = "store", dest = "sAccount", 
                    help = "Cascade project account,"
                           + " Ex: -accnt cascade12345")
parser.add_argument('-molcs', action = "store", dest = "sFileName1", 
                    help = "Molecule set with InchI codes,"
                           + " Ex: -molcs Molecule_Set.xlsx")
parser.add_argument('-metd', action = "store", dest = "sFileName2", 
                    help = "Mol, -metd Ex: Methods.xlsx")
parser.add_argument('-host', action = "store", dest = "sHost", 
                    help = "Cascade account name -host"
                           + " Ex: NetworkID@cascade.emsl.pnl.gov")
parser.add_argument('-dest', action = "store", dest = "sServerDir", 
                    help = "File directory in machine,"
                           + " Ex: -dest /dtemp/NetworkID/")
parser.add_argument('-nn', action = "store", dest = "iNodesNo", 
                    help = "Number of nodes, Ex: -nn 10", 
                    default = "1")
parser.add_argument('-nc', action = "store", dest = "iProcessorNoPerNode", 
                    help = "Number of processor per node, Ex: 16", 
                    default = "16")
parser.add_argument('-time', action = "store", dest = "sWalltime", 
                    help = "Maximum running time, Ex: -time 00:10:00", 
                    default = "00:05:00")

args = parser.parse_args()

#sLOCAL_DIR = args.sLocalDir
sACCOUNT = args.sAccount
sFILENAME1 = args.sFileName1
sFILENAME2 = args.sFileName2
sHOST = args.sHost
sSERVER_DIR = args.sServerDir
sNODES = args.iNodesNo
sPROCESSOR = args.iProcessorNoPerNode
sWalltime = args.sWalltime

sLOCAL_DIR = os.getcwd()

tools.CheckMustParser(sLOCAL_DIR, sFILENAME1, sFILENAME2, sHOST, sSERVER_DIR)

sLOCAL_DIR, sSERVER_DIR, \
    SHEET1, SHEET2 = exporters.ArrangeInputs(sLOCAL_DIR, sSERVER_DIR, 
                                             sFILENAME1, sFILENAME2)

# Store user id and password
scp.get_userinfo(sHOST)

tools.CheckInputMolecules(sLOCAL_DIR, SHEET1, SHEET2)

loggers.print_methods(SHEET2)
            
# Start simulation
#*************************************             
scp.simulate(SHEET1, SHEET2, sLOCAL_DIR, sACCOUNT, sHOST, sSERVER_DIR, 
             sNODES, sPROCESSOR, sWalltime)                 
#*************************************
#Simulation ends here 
            
exporters.isotropic_shieldings(SHEET1, SHEET2) 

exporters.chemical_shifts(SHEET1, SHEET2) 
      
loggers.error_file(SHEET1, SHEET2)

loggers.summary_file(SHEET1, SHEET2)

scp.clean_userinfo()  

end = time.time()
hours, rem = divmod(end-Start, 3600)
minutes, seconds = divmod(rem, 60)
print("\nIsicle Running Time: {:0>2}:{:0>2}:{:05.2f}".format(int(hours),
                                                             int(minutes),
                                                             seconds))  
