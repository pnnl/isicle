import subprocess
import glob
import logging
import numpy
import re
import csv
import pybel
import openbabel
import os


from timeit import default_timer as timer

logger = logging.getLogger(__name__)

class PkaMol:
    
     def __init__(self):
        self.acidic_atom_number=[]
        self.basic_atom_number=[]
        self.mol_id=[]
        self.min_acidic_pKa = []
        self.max_basic_pKa = []
        self.molfile=[]
        self.count = 0
        self.adduct_type=[]
        

    	     
     def writepka(self,fname,molids,adduct_types):
        with open(fname,'wb') as fp:
            a = csv.writer(fp,delimiter=',')
            a.writerow(["mol_id","molfile","acidic_atom_number",
                        "basic_atom_number","min_acidic_pKa","max_basic_pKa","adduct_type"])
                    
            for i in range(0,self.count):
                print self.mol_id[i]
                index = molids.index(self.mol_id[i])
                adduct_type = adduct_types[index]                
                
                out = [self.mol_id[i],self.molfile[i],str(self.acidic_atom_number[i]),
                   str(self.basic_atom_number[i]),str(self.min_acidic_pKa[i]),
                    str(self.max_basic_pKa[i]),str(adduct_type)]
                a.writerow(out)
     
     def readpka(self,fname):
         #fname = 'pka_values.csv'
         f = csv.reader(open(fname),delimiter = ',')
         line_num = 0
         for values in f:
             line_num += 1
             if(line_num ==1):
                 index0 = values.index('mol_id')
                 index1 = values.index('molfile')
                 index2 = values.index('acidic_atom_number')
                 index3 = values.index('basic_atom_number')
                 index4 = values.index('min_acidic_pKa')
                 index5 = values.index('max_basic_pKa')
                 index6 = values.index('adduct_type')
                
             if (line_num > 1):
                self.mol_id.append(values[index0])
                self.molfile.append(values[index1])
                self.acidic_atom_number.append(int(values[index2]))
                self.basic_atom_number.append(int(values[index3]))
                self.min_acidic_pKa.append(float(values[index4]))
                self.max_basic_pKa.append(float(values[index5]))
                self.adduct_type.append(values[index6])
        
         self.count= line_num - 1
         
        
     def cxcalcpka(self,pka_cmd,input_moldir,output_file,d_molids,d_adduct_types):
        """
        Uses Chemaxon cxcalc pKa plugin to determine the most acidic and 
        basic atoms in a molecule
        
        Runs on all '3D.mol' files present in the directory (input_moldir)
        
    
        """
        logger.info("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        logger.info("(PkaMol.cxcalcpka):cxcalc pKa calculations")
        logger.info("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
        counter=0
      #  fout = open(output_file, 'w')
        
        flist =input_moldir +  '/*_3D.mol'
        
        ftime = open('time_pkacalculations.csv','wb')
        ftime_write = csv.writer(ftime,delimiter=',')
        ftime_write.writerow(["mol_id","mol3D_file","time_seconds"])
            
        write_counter=0
        write_index = 0
        for mfile in glob.glob(flist):
           
    
            self.molfile.append(mfile)            
            counter = counter + 1
            self.count = counter
            write_counter = write_counter + 1
        
            s1 = re.split("/|,|\\\|\n",mfile)
            ns1 = len(s1)
            
            
            s2 = s1[ns1-1].split('_3D.mol')
            mol_id = s2[0]
            self.mol_id.append(mol_id)
            msg = "#"+str(counter)+": molecule id: " + mol_id
            logger.info(msg)
            
            
            msg = "mol file: "+ mfile
            logger.info(msg)
          
            #fout.write(msg+'\n')
         
            print pka_cmd
        
           # fout.write(cmd_line + '\n')
            cmd_line = pka_cmd + ',' + os.path.realpath(mfile)
            msg = "pKa command: " + cmd_line
            logger.info(msg)            
            
            
           
            args = cmd_line.split(",");
            print args
            
            ftime_start = timer()
            
            p=subprocess.Popen(args,stdout=subprocess.PIPE,shell=True)
            output =p.stdout.read()
            p_status = p.wait()
            print mfile
            print output
            a = output.split('\r\n')
            
            print a
            
            headers = a[0].split('\t')
            values =  a[1].split('\t')

            
            i_apKa1 = []
            i_apKa2 = []
            i_bpKa1 = []
            i_bpKa2 = []
            i_atoms= []    
            i_apKa1= headers.index("apKa1")
            i_apKa2 = headers.index("apKa2")
            i_bpKa1 = headers.index("bpKa1")
            i_bpKa2 = headers.index("bpKa2")
            i_atoms = headers.index("atoms")
            
            #for ih in headers:
      
       #         fout.write("%s," % ih)
                
        #    fout.write("\n")
            
         #   for ih in values:
                
          #      fout.write("%s,"% ih)
                
           # fout.write("\n")
           # fout.write("----------------------------------------------\n")
            
                 
            print len(values)     
            if(values[i_apKa1]!="pka:FAILED") and len(values) > 2:
                
                atom_nums = values[i_atoms].split(',')
                msg= ["atom numbers: ", atom_nums]
                logger.info(msg)
            
            #print "%s, %s, %s, %s, %s" %(i_apKa1,i_apKa2,i_bpKa1,i_bpKa2,i_atoms)
            
           
                nacidpka = 0
                nbasicpka = 0
                index_acidic_atom = []
                index_basic_atom = []
                vs = filter(None,values[i_apKa1:i_apKa2+1])
                if vs: # list is not empty
                    msg = ["Acidic pKas: ",vs]
                    logger.info(msg)
                    
                    nacidpka =len(vs) # get number of acidic pKa values
                    num_vs = numpy.array(vs,dtype=float) # convert str to float
                    
                    # select only those that have a hydrogen atom attached 
                    
                    min_num_vs = numpy.min(num_vs)
                    self.min_acidic_pKa.append(min_num_vs) # lowest value
                    # get index corresponding to the lowest acidic pKa in 
                    index_acidic_atom = numpy.argmin(num_vs) # get the index of the minimum value     
                    alternate_index_acidic_atom = numpy.argmax(num_vs)
                    msg = "Lowest acidic pKa = " + str(min_num_vs)
                    logger.info(msg)
                    
                    
                    atoms_acidpka = atom_nums[0:nacidpka]
                    i1 = numpy.int(atoms_acidpka[index_acidic_atom])
                    
                    mol = pybel.readfile("mol",mfile).next()
          
                    iatom = mol.atoms[i1-1]
                    nbatoms = openbabel.OBAtomAtomIter(iatom.OBAtom)
                    
                    found = 0
                    for nb in nbatoms:
                        nb_atomnum = nb.GetAtomicNum()
                        if(nb_atomnum==1):
                            found = 1
                    if (found == 0):
                        msg = "Atom number with lowest acidic pKa has no hydrogen atom attached"
                        logger.info(msg)
                        msg = "Checking the other site"
                        logger.info(msg)
                        
                        i1 = numpy.int(atoms_acidpka[alternate_index_acidic_atom])
                        iatom = mol.atoms[i1-1]
                        nbatoms = openbabel.OBAtomAtomIter(iatom.OBAtom)
                        
                        for nb in nbatoms:
                            nb_atomnum = nb.GetAtomicNum()
                            if(nb_atomnum==1):
                                found = 1
                    if (found == 1):
                        
                        self.acidic_atom_number.append(i1)
                        msg = "Atom number selected for deprotonation = "
                        msg = msg + numpy.str(i1) + '\n'
                        logger.info(msg)
                    
                    else: 
                        msg = "No acidic site found"
                        logger.info(msg)
                        self.min_acidic_pKa.append(10000) # a dummy value
                        self.acidic_atom_number.append(-1)
                else:
                    self.min_acidic_pKa.append(10000) # a dummy value
                    self.acidic_atom_number.append(-1)
                    
                vs = filter(None,values[i_bpKa1:i_bpKa2+1])
                if vs:  # list is not empty
                    msg = ["Basic pKas: ",vs]
                    logger.info(msg)
                    nbasicpka = len(vs) # get number of basic pKa values
                    num_vs = numpy.array(vs,dtype=float)
                    max_num_vs = numpy.max(num_vs)
                    self.max_basic_pKa.append(max_num_vs) # highest value
                    index_basic_atom = numpy.argmax(num_vs) # get the index of the maximum value
                    msg = "Highest basic pKa:" + str(max_num_vs)
                    logger.info(msg)
                    
                    npka = nacidpka + nbasicpka
                    
                    atoms_basicpka = atom_nums[nacidpka:npka]
                    i1 = numpy.int(atoms_basicpka[index_basic_atom])
                    self.basic_atom_number.append(i1)
                    
                    msg = "Atom number with the highest basic pKa = "
                    msg = msg + numpy.str(i1) + '\n'
                    logger.info(msg)
                else:
                    self.max_basic_pKa.append(10000)
                    self.basic_atom_number.append(-1)
            else:
                msg = "pKa calculation failed."
                logger.info(msg)
                self.basic_atom_number.append(-1)
                self.max_basic_pKa.append(10000)
                self.min_acidic_pKa.append(10000) # a dummy value
                self.acidic_atom_number.append(-1)
            ftime_end = timer()
                                
           
            out = [mol_id,mfile,str(ftime_end - ftime_start)]
            ftime_write.writerow(out)
            
            if (write_counter == 100):
                write_counter = 0
                write_index = write_index + 1
                fout = 'pka_values' + numpy.str(write_index) +'.csv'
                self.writepka(fout,d_molids,d_adduct_types)
        ftime.close()    
       # fout.close()
        msg = '(PkaMol.cxcalcpka): cxcalc pKa calcualations completed.'
        logger.info(msg)

############################################################################
# Function to read InChI list
