from argparse import ArgumentParser, FileType
from pandas import DataFrame, read_table
from chembl_ikey import inchi_to_inchikey
from os import path
from pathlib import Path
import yaml

# Run file from command line with location of inchi source files, optional naming
#   ex: python generate_input.py ./processed_InChI.txt
#           default file output ./input/inchikey.inchi
#   or: python generate_input.py ./processed_InChI.txt --unique --output
#           to write unique filepaths and a unique output path


def main():
    if args.output is not None:
        config = args.output
    else:
        stream = open('config.yaml', 'r')
        configuration = yaml.load(stream)
        config = path.join(configuration['path'], './input/')
        stream.close()
    unique_inchis() if args.unique else default_inchis(config)


def default_inchis(config):
    i = 0
    for row in df.values:
        key = inchi_to_inchikey(row[0])
        filename = path.normpath(path.join(config, '%s.inchi' % key))
        checkingfile(i, filename)
        i += 1


def unique_inchis():
    i = 0
    for row in df.values:
        key = inchi_to_inchikey(row[0])
        filename = os.path.normpath(input("Type filepath including .inchi for %s \n" % key))
        checkingfile(i, filename)
        i += 1


def checkingfile(i, filename):
    try:
        abs_path = Path(filename).resolve(strict=True)
    except FileNotFoundError:
        df.iloc[i].to_csv(filename, index=False, sep='\t', header=None)
    except RuntimeError:
        print("Invalid. File not saved.\n")
    else:
        overwrite = input("%s already exists. Overwrite <y/n>? \n" % filename)
        if overwrite == 'y':
            df.iloc[i].to_csv(filename, index=False, sep="\t", header=None)
        elif overwrite == 'n':
            new_filename = os.path.normpath(input("Type new filepath, including .inchi\n"))
            checkingfile(i, new_filename)
        else:
            print("Invalid. File not saved.\n")


if __name__ == '__main__':
    parser = ArgumentParser(description="A program to convert a file of inchis to individual inchikey.inchi files.")
    parser.add_argument('filepath', type=FileType('r'), help="Include valid filepath to inchi file.")
    parser.add_argument('--unique', nargs='?', const='store_true', help="Include unique to write uniquefilename.inchi")
    parser.add_argument('--output', nargs='?', const=None,
                        help="Allows for unique inchi output location (Default path: config.yaml[''path'']/input/")
    args = parser.parse_args()
    df = read_table(args.filepath, sep='\n', header=None)
    df = df.replace('"', '')
    main()
    parser.exit()
