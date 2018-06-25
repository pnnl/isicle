from argparse import ArgumentParser,FileType
from pandas import DataFrame, read_table
from chembl_ikey import inchi_to_inchikey
###chembl_ikey opensource github
### file ikey.py slightly modified, added to line 27: bStdFormat = 0
import os

### Run file from command line with location of inchi source files
#   ex: python generate_input.py ./processed_InChI.txt

parser = ArgumentParser()
parser.add_argument('filepath',type=FileType('r'))
args = parser.parse_args()
print(args.filepath)
df = read_table(args.filepath,sep='\n',header=None)
df = df.replace('"','')

def main():
    naming = input("Write unique filepath for .inchi file? Default ./input/InChikey.inchi <y/n> \n")
    if naming == 'y': unique_inchis()
    else: default_inchis()
    parser.exit()

def unique_inchis():
    i=0
    for row in df.values:
        key = inchi_to_inchikey(row[0])
        filename = os.path.normpath(input("Write filepath including .inchi for %s \n" % key))
        file_present = os.path.isfile(filename)
        if file_present:
            print(OSError)
            overwrite = input("%s already exists. Overwrite <y/n>? \n" % filename)
            if overwrite == 'y':
                df.iloc[i].to_csv(filename,index=False,sep="\t")
            elif overwrite =='n':
                new_filename = os.path.normpath(input('Type new filename.\n'))
                df.iloc[i].to_csv(new_filename,index=False,sep="\t")
            else:
                print("Invalid. File not saved.\n")
        else:
            df.iloc[i].to_csv(filename,index=False,sep="\t")
        i+=1
def default_inchis():
    i=0
    for row in df.values:
        key = inchi_to_inchikey(row[0])
        filename = "./input/%s.inchi" % key
        file_present = os.path.isfile(filename)
        if file_present:
            print(OSError)
            overwrite = input("%s already exists. Overwrite <y/n>? \n" % filename)
            if overwrite == 'y':
                df.iloc[i].to_csv(filename,index=False,sep="\t")
            elif overwrite =='n':
                new_filename = os.path.normpath("./input/%s.inchi" % input('Type new filename. \n'))
                df.iloc[i].to_csv(new_filename,index=False,sep="\t")
            else:
                print("Invalid. File not saved.\n")
        else:
            df.iloc[i].to_csv(filename,index=False,sep="\t")
        i+=1

if __name__ == '__main__':
    main()
