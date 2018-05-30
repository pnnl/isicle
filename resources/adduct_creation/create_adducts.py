# -*- coding: utf-8 -*-
"""
Created on Sat Sep 03 01:00:25 2016

@author: Dennis Thomas
"""

import sys
sys.path.insert(0, '../')
sys.path.insert(0,'../../')

import geometry
import pka
import readio
import csv
import pybel
import openbabel
import datetime
import logging
import glob
import numpy


from timeit import default_timer as timer



def main():
    # read inchi file list
    
    # read input parameters from file, parameters.in
    
    dt = datetime.datetime.now()
    dtstr=dt.strftime("%A, %B %d, %Y %I:%M %p")
    
    # create a name for the log file based on start date and time
    
    logfile='run_' + str(dt.month) + "_" + str(dt.day) + "_" + str(dt.year)+ "_"
    logfile = logfile + str(dt.hour) + "_" + str(dt.minute) + "_" + str(dt.second) + ".log"
    print logfile
    
    logging.basicConfig(filename=logfile,level=logging.INFO)
    logger = logging.getLogger(__name__)
    msg = "Started " + dtstr
    logger.info(msg)
    
    fname = 'parameters.in'  
    msg = "Reading parameters from " + fname    
    logger.info(msg)
    
    uio = readio.IOFile(fname)
    uio.readFile()
    
#    # read inCHi list
#    msg = "Reading InChI list"
#    logger.info(msg)
#    d= geometry.InChI(uio.inchilist_dir,uio.inchilist_file)
#    d.read_inChiFile()
#    
#    #Convert inChI to .mol format
#    time_inchi2molxyz_seconds = []
#    mol3Dfiles = []
#    for i in range(0,d.num_cmpnds):
#        tstart = timer()
#       
#        inchi_str = d.inchi[i]
#        cpd_id = d.object_id[i]
#       
#        (mol,molfile,mol2Dfile,mol3Dfile) = geometry.inchi2mol(inchi_str,cpd_id,
#                                        "/",uio.inchi2mol_dir,uio.forcefield)
#    
#    # Create .xyz. file
#        geometry.createXYZFile(mol,mol3Dfile,uio.inchi2xyz_dir)
#        tend = timer()
#        mol3Dfiles.append(mol3Dfile)
#        time_inchi2molxyz_seconds.append(tend-tstart)
#        
#    with open('time_inch2molxyz.csv','wb') as fp:
#            a = csv.writer(fp,delimiter=',')
#            a.writerow(["mol_id","mol3D_file","time_seconds"])
#                    
#            for i in range(0,d.num_cmpnds):
#                out = [d.object_id[i],mol3Dfiles[i],str(time_inchi2molxyz_seconds[i])]
#                a.writerow(out)
#            
#    fp.close()

# calculate pKa (acidity increases with decrease in pKa (i.e., more negative))
    #msg = "Starting pKa calculations"
    #logger.info(msg)
    
    
    pkamol = pka.PkaMol()
    #pkamol.cxcalcpka(uio.pka_cmd_line,uio.inchi2mol_dir,uio.pka_outputfile)
    #pkamol.writepka('pka_values.csv')


    logger.info("Creating adduct files")

    basic_site_adduct_types=["(+H)","(+Na)"]
    basic_site_adduct_elements = ["H","Na"]
    
    basic_site_adduct_bond = ["single","none"]
    basic_site_adduct_bonded = ["yes","no"]
    acidic_site_adduct_types=["(+De)"]
    
    
    flist = 'pka_values*'
    fcount = 0    
    for mfile in glob.glob(flist):
        msg = "Reading pKa values"    
        logger.info(msg)    
        pkamol.readpka(mfile)
    
        print "count = " , pkamol.count
        fcount = fcount + 1
        f1 = 'time_adductcreation' + numpy.str(fcount) + '.csv'
        ftime = open(f1,'wb')
        ftime_write = csv.writer(ftime,delimiter=',')
        ftime_write.writerow(["mol_id","mol3D_file","adduct_createtime_sec"])
 #       ftime_write.writerow(["mol_id","mol3D_file","protonation_time_sec",
  #                            "sodiation_time_sec","deprotontation_time_sec"])
    # add hydrogen (to the most basic atom)
        for i in range(0,pkamol.count): # for each molecule, create the adducts
                   
            msg = "For molecule id " + pkamol.mol_id[i]
            logger.info(msg)
            msg = ".mol file: " + pkamol.molfile[i]
            logger.info(msg)
            
            
            # create protonated adduct
         #   time_prot = numpy.nan
         #   time_sod = numpy.nan
          #  time_deprot = numpy.nan
            time_adduct = numpy.nan
            if(pkamol.adduct_type[i] in basic_site_adduct_types):
                index = basic_site_adduct_types.index(pkamol.adduct_type[i])
                element=basic_site_adduct_elements[index]
                bond_type=basic_site_adduct_bond[index]
                bonded = basic_site_adduct_bonded[index]
    
                if(pkamol.basic_atom_number[i]!= -1):
                    adduct = []
                    xyzr = []
                    tstart = timer()            
                    mol = pybel.readfile("mol",'../' + pkamol.molfile[i]).next()
                    print pkamol.molfile[i]
                    print str(len(mol.atoms))
       # (adduct,xyzr) = geometry.addAtomToMol(mol,'H',
        #                pkamol.basic_atom_number[i],"single","yes")
                    (adduct,xyzr) = geometry.addAtomToMol(mol,element,
                                    pkamol.basic_atom_number[i],bond_type,bonded)              
                    adduct_molfile = uio.adductmol_dir + '/' + pkamol.mol_id[i] + "+" + element + '.mol'
                    adduct_xyzfile = uio.adductxyz_dir + '/' + pkamol.mol_id[i] + "+" + element + '.xyz'
                    adduct_pngfile = uio.adductmol_dir + '/' + pkamol.mol_id[i] + "+" + element + '.png'
                    adduct_mol2file = uio.adductmoltwo_dir + '/' + pkamol.mol_id[i] + "+" + element + '.mol2'   
	    # optimization
                    print "Optimizing the protonated molecule..."
                    adduct.localopt(forcefield=uio.forcefield,steps=1000)
        
                    adduct.draw(show=False,filename=adduct_pngfile,usecoords=False,update=False)
        
        
                    adduct.write("xyz",adduct_xyzfile,True)
                    tend = timer()
                    time_adduct = tend - tstart
       # adduct= pybel.readfile("xyz",adduct_xyzfile).next()
#                    adduct.calccharges()                    
#                    adduct.write("mol",adduct_molfile,True)
#                    adduct.write("mol2",adduct_mol2file,True)
        #print "Finished writing out .mol, .mol2, and .xyz files"            
                    msg = "Created adduct molecule by adding " + element + " at " 
                    msg = msg + str(pkamol.basic_atom_number[i])
                    logger.info(msg)
        
                else:
                    msg = "No basic site found to create protonated (+H) and sodiated (+Na) adducts"
                    logger.info(msg)


    # create de-protonated adduct
            #print "adduct type: " + pkamol.adduct_type[i]
            if(pkamol.adduct_type[i] in acidic_site_adduct_types):
                #print "creating adduct" + pkamol.adduct_type[i]
                index = acidic_site_adduct_types.index(pkamol.adduct_type[i])
                if(pkamol.acidic_atom_number[i]!= -1):
    
                    adduct = []
                    tstart = timer()
                    mol = pybel.readfile("mol",'../' + pkamol.molfile[i]).next()
  
                    iatom = mol.atoms[pkamol.acidic_atom_number[i]-1]
    # get the neighboring atoms of the selected atom
                    nbatoms = openbabel.OBAtomAtomIter(iatom.OBAtom)
    # get the atom numbers of the neighboring atoms
                    h_atom_index = []
                    found = 0
    
                    for nb in nbatoms:
                        nb_atomnum = nb.GetAtomicNum()
                        if(nb_atomnum==1 and found == 0): # if the neighboring atom is hydrogen
                            h_atom_index.append(nb.GetId())
                            msg = "acidic atom index = " + str(pkamol.acidic_atom_number[i]-1)
                            logger.info(msg)
                            msg = "hydrogen atom index to remove" + str(h_atom_index[0])
                            logger.info(msg)
                            print h_atom_index
                            found = 1
                    if(found == 1):
                        (adduct,total_chg)=geometry.removeAtomFromMol(mol,h_atom_index[0])
 
                        adduct_molfile = uio.adductmol_dir + '/'+ pkamol.mol_id[i] + '+De.mol'
                        adduct_xyzfile = uio.adductxyz_dir + '/' + pkamol.mol_id[i] + '+De.xyz'
                        adduct_pngfile = uio.adductmol_dir + '/' + pkamol.mol_id[i] + '+De.png'        
                        adduct_mol2file = uio.adductmoltwo_dir + '/'+ pkamol.mol_id[i] + '+De.mol2'

                        print "Optimizing the de-protonated molecule..."
                        adduct.localopt(forcefield=uio.forcefield,steps=1000)
                        adduct.draw(show=False,filename=adduct_pngfile,usecoords=False,update=False)
                        adduct.write("xyz",adduct_xyzfile,True)
    #    adduct= pybel.readfile("xyz",adduct_xyzfile).next()
#                        adduct.calccharges()                        
#                        adduct.write("mol",adduct_molfile,True)
#                        adduct.write("mol2",adduct_mol2file,True)

                        tend  = timer()
                        time_adduct = tend - tstart
                        msg = "De-protonated molecule by removing hydrogen "
                        msg = msg + " (atom index = " + str(h_atom_index[0]) + "),\n"
                        msg = msg + "attached to atom number" 
                        msg = msg + str(pkamol.basic_atom_number[i])
                        logger.info(msg)
                else:
                    msg = "No acidic site found to create de-protonated adducts"
                    logger.info(msg)       
            out = [pkamol.mol_id[i],pkamol.molfile[i],str(time_adduct)]
            ftime_write.writerow(out)          
   # out = [pkamol.mol_id[i],pkamol.molfile[i],str(time_prot),
    #       str(time_sod),str(time_deprot)]
        
        ftime.close()
                
               
            
#==============================================================================
#             
#             if(found==0):
#                 msg = "No hydrogen atom available at selected basic site"
#                 logger.info(msg)        
#         else:
#             msg = "No acidic site found to create de-protonated adducts"
#             logger.info(msg)
#==============================================================================
    
        
if __name__ == '__main__':
    main()
            
        
