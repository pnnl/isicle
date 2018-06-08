import logging
import os

logger = logging.getLogger(__name__)

class IOFile:
    
    """
    This class stores the user-specified input/output file and directory names .
    """
    def __init__(self,inputfile):
        self.input_file = inputfile # input file
        self.inchilist_dir = []
        self.inchilist_file = []
        self.inchi2mol_dir = []
        self.inchi2xyz_dir = []
        self.adductxyz_dir = []
        self.adductmol_dir = []
        self.adductmoltwo_dir = []
        self.forcefield = []
        self.pka_cmd_line = []
        self.pka_outputfile = []
    
    def readFile(self):
   
        f = open(self.input_file, 'r')
        print "input parameter file:",self.input_file
        
        for line in f:
            
            values = line.split(",",1)
           
            if len(values) == 2:
                if(values[0]=="INCHILIST_FILEDIR"):
                    self.inchilist_dir = values[1].strip()
                    logger.info("(IoFile.readFile): InChI list read from file in directory: " 
                    + self.inchilist_dir)
                elif(values[0]=="INCHILIST_FILENAME"):
                    self.inchilist_file = values[1].strip()
                    logger.info("(IoFile.readFile): InChI list read from file: " + self.inchilist_file 
                    + '\n')
                elif(values[0]=="INCHI2MOL_DIR"):
                    self.inchi2mol_dir = values[1].strip()
                    logger.info("(IoFile.readFile): .mol files derived from InChIs will be saved in directory : " 
                    + self.inchi2mol_dir)
                    if os.path.exists(self.inchi2mol_dir)==0:
                        os.mkdir(self.inchi2mol_dir)
                elif(values[0]=="INCHI2XYZ_DIR"):
                    self.inchi2xyz_dir = values[1].strip()
                    logger.info("(IoFile.readFile): .xyz files derived from InChIs will be saved in directory : " 
                    + self.inchi2xyz_dir)
                    if os.path.exists(self.inchi2xyz_dir)==0:
                        os.mkdir(self.inchi2xyz_dir)
                        
                elif(values[0]=="ADDUCTXYZ_DIR"):
                    self.adductxyz_dir = values[1].strip()
                    logger.info("(IoFile.readFile): Adduct .xyz files will be saved in directory: " 
                    + self.adductxyz_dir)
                    if os.path.exists(self.adductxyz_dir)==0:
                        os.mkdir(self.adductxyz_dir)
                        
                elif(values[0]=="ADDUCTMOL_DIR"):
                    self.adductmol_dir = values[1].strip()
                    logger.info("(IoFile.readFile): Adduct .mol files will be saved in directory: "
                    + self.adductmol_dir)
                    if os.path.exists(self.adductmol_dir)==0:
                        os.mkdir(self.adductmol_dir)
                
                elif(values[0]=="ADDUCTMOLTWO_DIR"):
                    
                    self.adductmoltwo_dir = values[1].strip()
                    logger.info("(IoFile.readFile): Adduct .mol2 files will be saved in directory: "
                    + self.adductmoltwo_dir)
                    if os.path.exists(self.adductmoltwo_dir)==0:
                        os.mkdir(self.adductmoltwo_dir)

                elif(values[0]=="FF_GEOMETRY_OPTIMIZATION"):
                    self.forcefield = values[1].strip()
                    logger.info("(IoFile.readFile): Forcefield used for 2D geometry optimization: " 
                    + self.forcefield)
                            
                elif(values[0]=="PKA_COMMAND_LINE"):
                    self.pka_cmd_line = values[1].strip()
                    logger.info("(IoFile.readFile): The command used for pKa calculations: " 
                    + self.pka_cmd_line)
                
                elif(values[0]=="PKA_OUTPUT_FILE"):
                    self.pka_outputfile = values[1].strip()
                    logger.info("(IoFile.readFile): pKa values will be saved in file: " 
                    + self.pka_outputfile)
                    
        f.close()
        
