from isicle.interfaces import MolecularStringInterface
from rdkit import Chem
from rdkit.Chem import SaltRemover, AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
import subprocess
import argparse
import sys

"""
Methods:

load(self,path):
    Loads data from filepath.
desalt(self):
    Converts SMILES to RDKit mol object.
    Desalts RDKit mol object using Chem.SaltRemover module.
    Returns desalted RDKit SMILES string.
    (not instituted) Salts to be removed with a config file.
neutralize(self):
    Converts SMILES to RDKit mol object.
    Neutralizes SMILES string.
    Returns neutralized SMILES string.
tautomerize(self,alt=False):
    Converts SMILES to RDKIT mol object.
    Generate tautomers according to RDKit TautomerEnumerator() method.
    Default returns first tautomer generated.
    If alt=True, returns all generated tautomers
to_smiles(self,type='mol'):
    Default converts RDKit mol object to RDKit canonical SMILES.
    If type='SMILES', converts SMILES to RDKit canonical SMILES.
to_mol(self,type='SMILES'):
    Default converts SMILES string to RDKit mol object.
    If type='SMARTS', converts SMARTS string to RDKit mol object.
"""


class MolecularString(MolecularStringInterface):
    """Desalt, neutralize, and tautomerize SMILES"""

    def __init__(self):
        self.contents = None

    def load(self, path: str):
        # offer commandline string, not necessarily path
        with open(path, 'r') as f:
            self.contents = f.readlines()
        return self.contents

    def desalt(self):
        """Generate desalted RDKit SMILES from SMILES string"""
        remover = SaltRemover.SaltRemover(defnFormat='smiles', defnData="[Cl]")
        # defnData="[Cl,Br,Na]" *sample definition of salts to be removed*
        # add iterator for salts listed in config?
        # set potential salts to be removed in a config file

        try:
            res, deleted = remover.StripMolWithDeleted(self.to_mol(self.contents))
            res = self.to_smiles(res)
            # using StripMolWithDeleted instead of StripMol
            # add functionality to track removed salts
            # atomno = res.GetNumAtoms
            # if relevant to future use, returns atom count post desalting
            return res

        except NotImplementedError

    def neutralize(self):
        """Generate neutralized RDKit SMILES from SMILES string"""

        def _InitializeNeutralisationReactions():
            patts = (
                # Imidazoles
                ('[n+;H]', 'n'),
                # Amines
                ('[N+;!H0]', 'N'),
                # Carboxylic acids and alcohols
                ('[$([O-]);!$([O-][#7])]', 'O'),
                # Thiols
                ('[S-;X1]', 'S'),
                # Sulfonamides
                ('[$([N-;X2]S(=O)=O)]', 'N'),
                # Enamines
                ('[$([N-;X2][C,N]=C)]', 'N'),
                # Tetrazoles
                ('[n-]', '[nH]'),
                # Sulfoxides
                ('[$([S-]=O)]', 'S'),
                # Amides
                ('[$([N-]C=O)]', 'N'),
            )
            return [(self.to_mol(x, type='SMARTS'), self.to_mol(y)) for x, y in patts]
        try:
            reactions = _InitializeNeutralisationReactions()

            mol = self.to_mol(self.contents)
            replaced = False
            for i, (reactant, product) in enumerate(reactions):
                while mol.HasSubstructMatch(reactant):
                    replaced = True
                    rms = AllChem.ReplaceSubstructs(mol, reactant, product)
                    mol = rms[0]
            if replaced:
                return self.to_smiles(mol)
            else:
                return self.to_smiles(self.contents, type='SMILES')

        except:
            raise NotImplementedError

    def tautomerize(self, alt=False):
        """Generate tautomerized RDKit mol object(s) from SMILES string"""
        # source: https://rdkit.blogspot.com/2020/01/trying-out-new-tautomer.html
        # Discuss noted double bond changes
        try:
            enumerator = rdMolStandardize.TautomerEnumerator()

            mol = self.to_mol(self.contents)
            res = [self.contents]
            tauts = enumerator.Enumerate(mol)
            smis = [self.to_smiles(x) for x in tauts]
            s_smis = sorted((x, y)
                            for x, y in zip(smis, tauts) if x != self.contents)
            res += [y for x, y in s_smis]

            if alt is True:
                return res
            else:
                return res[0]
        except:
            raise NotImplementedError

    def to_smiles(self, type='mol'):
        """Convert RDKit mol object or SMILES to RDKit canonical SMILES"""
        try:
            if type == 'mol':
                res = MolToSmiles(self.contents)
                return res
            elif type == 'SMILES':
                res = MolToSmiles(MolFromSmiles(self.contents))
        except:
            raise NotImplementedError

    def to_mol(self, type='SMILES'):
        """Convert SMILES or SMARTS string to RDKit mol object"""
        try:
            if type == 'SMILES':
                res = MolFromSmiles(self.contents)
                return res
            elif type == 'SMARTS':
                res = MolFromSmarts(self.contents)
                return res
        except:
            raise NotImplementedError

    def to_file(self, type, fn):
        """Write contents to file"""
        if type == 'string':
            with open(path, 'w') as f:
                f.write(self.contents)
                f.close()
        elif type == 'pdb':
            Chem.rdmolfiles.MolToPDBFile(self.contents, fn)
