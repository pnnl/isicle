from isicle.interfaces import GeometryGenerationInterface
from rdkit import Chem
from rdkit.Chem import AllChem
import pybel
import pickle
import numpy as np

# TO DO, read from xyz or mol2 to yield class
# in utils.py, read_mol, Mol, pop_aom, push_atom functions


class GeometryGeneration(GeometryGenerationInterface):
    def __init__(self):
        self.contents = None

    def load(self, path: str):
        """Load in the data file"""
        with open(path, 'r') as f:
            self.contents = f.readlines()
        self.path = path
        return self.contents

    def inputto2D(self, path: str):
        """Convert SMILES/INCHI into 2D geometry using RDKit"""
        if self.path.endswith('.inchi'):
            mol = Chem.MolFromInchi(self.contents)

        if self.path.endswith('.smi'):
            mol = Chem.MolFromSmiles(self.contents)

        molH = Chem.AddHs(mol)
        self.molecule2D = Chem.MolToMolBlock(molH)

        if self.path.endswith('.xyz'):
            """Convert XYZ to Mol"""
            xyz = next(pybel.readfile("xyz", FILE))
            name = ((self.path).split('.')[0]).split('/')[-1] + '.mol'
            mol = xyz.write('mol', name, erwrite=True)

        if self.path.endswith('.mol') or self.path.endswith('.mol2'):
            self.molecule2D = self.contents

        file1 = open('geom_2D.mol', 'w+')
        file1.write(self.molecule2D)
        file1.close()
        return self.molecule2D

    def convert3D(self, path: str):
        """Convert 2D geometry mol file into 3D geometry using RDKit"""

        self.molecule3D = AllChem.MMFFOptimizeMolecule(self.molecule2D)
        file2 = open('geom_3D.mol', 'w+')
        file2.write(self.molecule3D)
        file2.close()
        return self.molecule3D

    def total_partial_charge(self):
        self.charge = np.array([a.partialcharge for a in self.atoms]).sum()
        return self.charge

    def natoms(self):
        self.number = len(self.atoms)
        return self.number

    def save(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        return
