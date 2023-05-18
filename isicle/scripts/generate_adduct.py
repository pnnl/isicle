import argparse
import os
from openbabel import openbabel as ob
from os.path import abspath, join
from isicle.resources import geometry
from isicle.utils import read_pka, read_mol, write_string
from isicle import __version__

# import subprocess
# from os.path import *


def parse_mol_index(mol):
    hBondList = []
    nonHList = []
    for obatom in ob.OBMolAtomIter(mol.OBMol):
        # for obatom in mol.atoms:
        if obatom.GetAtomicNum() not in [1, 6]:  #!= 1:
            idx = obatom.GetIdx()
            nonHList.append(idx)
            # nbatoms = ob.OBAtomAtomIter(obatom.OBAtom)
            for ob2atom in ob.OBAtomAtomIter(obatom):
                if ob2atom.GetAtomicNum() == 1:
                    hBondList.append(idx)
    return (nonHList, list(set(hBondList)))


def create_adduct(mol, adduct, idx, forcefield="mmff94", steps=500):
    if adduct == "-H":
        hidx = geometry.nearestHydrogen(mol, idx)
        mol = geometry.removeAtomFromMol(mol, hidx)
    elif "+" in adduct:
        atom = adduct.split("+")[-1]
        if atom == "Na":
            mol = geometry.addAtomToMol(mol, atom, idx, covalent=False)
        elif atom == "H":
            mol = geometry.addAtomToMol(mol, atom, idx, covalent=True)

    # _builder = ob.OBBuilder()
    # _builder.Build(mol.OBMol)
    mol.localopt(forcefield=forcefield, steps=steps)

    # # adjust partial charge
    # if adduct == '-H':
    #     mol.atoms[idx].OBAtom.SetPartialCharge(mol.atoms[idx].partialcharge - 1)

    return mol


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate adduct from base geometry")
    parser.add_argument("infile", help="path to input .mol2 file")
    parser.add_argument("adduct", help="adduct type")
    parser.add_argument(
        "-ff", "--forcefield", type=str, default="gaff", help="forcefield type"
    )
    parser.add_argument(
        "-s",
        "--steps",
        type=int,
        default=500,
        help="number of forcefield optimization steps",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=__version__,
        help="print version and exit",
    )

    args = parser.parse_args()
    args.adduct = args.adduct[1:-1]

    # read input
    mol = read_mol(args.infile, fmt="mol2")

    if args.adduct == "neutral":
        molList = [mol]
        addSiteList = [0]
        pass
    else:
        nonHList, hBondList = parse_mol_index(mol)
        molList = []
        # generate adduct
        if "+" in args.adduct:
            addSiteList = nonHList
            for atomSite in nonHList:
                # reset mol object
                mol = read_mol(args.infile, fmt="mol2")
                molList.append(
                    create_adduct(
                        mol,
                        args.adduct,
                        atomSite,
                        forcefield=args.forcefield,
                        steps=args.steps,
                    )
                )
        elif "-" in args.adduct:
            addSiteList = hBondList
            for atomSite in hBondList:
                # reset mol object
                mol = read_mol(args.infile, fmt="mol2")
                molList.append(
                    create_adduct(
                        mol,
                        args.adduct,
                        atomSite,
                        forcefield=args.forcefield,
                        steps=args.steps,
                    )
                )

    for molObj, atomSite in zip(molList, addSiteList):
        molID = args.infile.split("/")[-1].split(".")[0]
        fp = abspath(join("output", "adducts", f"geometry_{args.adduct}", f"{molID}"))
        if not os.path.exists(fp):
            os.makedirs(fp)
        fp = join(fp, f"{atomSite}")
        molObj.write("mol2", fp + ".mol2", overwrite=True)
        molObj.write("xyz", fp + ".xyz", overwrite=True)
        molObj.write("pdb", fp + ".pdb", overwrite=True)

        # reassign charge
        if args.adduct in ["+H", "+Na"]:
            # tmp = splitext(args.mol2)[0] + '.tmp.mol2'
            # subprocess.call('obabel %s -O %s --partialcharge eem' % (args.mol2, tmp), shell=True)
            # subprocess.call('mv %s %s' % (tmp, args.mol2), shell=True)
            write_string("1.0", fp + ".charge")

        elif args.adduct == "-H":
            write_string("-1.0", fp + ".charge")

        elif args.adduct == "neutral":
            write_string("0.0", fp + ".charge")
