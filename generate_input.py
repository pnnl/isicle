from argparse import ArgumentParser,FileType
from pandas import DataFrame, read_table
from chembl_ikey import inchi_to_inchikey
###chembl_ikey opensource github
### file ikey.py slightly modified, added to line 27: bStdFormat = 0
import os
from pathlib import Path

### Run file from command line with location of inchi source files, optional naming
#   ex: python generate_input.py ./processed_InChI.txt
#           default file output ./input/inchikey.inchi
#   or: python generate_input.py ./processed_InChI.txt --unique
#           to write unique filepaths

parser = ArgumentParser()
parser.add_argument('filepath',type=FileType('r'))
parser.add_argument('--unique',nargs='?')
args = parser.parse_args()
df = read_table(args.filepath,sep='\n',header=None)
df = df.replace('"','')

def main():
    unique_inchis() if args.unique else default_inchis()

def default_inchis():
    i=0
    for row in df.values:
        key = inchi_to_inchikey(row[0])
        filename = "./input/%s.inchi"%key
        checkingfile(i,filename)
        i+=1

def unique_inchis():
    i=0
    for row in df.values:
        key = inchi_to_inchikey(row[0])
        filename = os.path.normpath(input("Write filepath including .inchi for %s \n" % key))
        checkingfile(i,filename)
        i+=1

def checkingfile(i,filename):
    try:
        abs_path = Path(filename).resolve(strict=True)
    except FileNotFoundError:
        df.iloc[i].to_csv(filename,index=False,sep='\t',header=None)
    except RuntimeError:
        print("Invalid. File not saved.\n")
    else:
        overwrite = input("%s already exists. Overwrite <y/n>? \n" % filename)
        if overwrite == 'y':
            df.iloc[i].to_csv(filename,index=False,sep="\t",header=None)
        elif overwrite =='n':
            new_filename = os.path.normpath(input('Type new filename.\n'))
            checkingfile(i,new_filename)
        else:
            print("Invalid. File not saved.\n")

if __name__ == '__main__':
    main()
    
parser.exit()
