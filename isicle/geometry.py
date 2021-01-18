from isicle.interfaces import GeometryGenerationInterface
from rdkit import Chem
from rdkit.Chem import AllChem
import pybel
import pickle
import numpy as np

# TO DO, read from xyz or mol2 to yield class
# in utils.py, read_mol, Mol, pop_aom, push_atom functions


def load_pickle():
    # return autodetect (isinstance)


def load_xyz():
    # return Geometry instance


def load_mol()
    # return Geometry instance


def load_mol2()
    # return Geometry instance


class Geometry(GeometryInterface):

    def to_2D(self):
        # return Geometry2D instance

    def to_3D(self):
        # return Geomtery3D instance

    def save(self):
        # pickle save


class Geometry2D(Geometry):
    '''
    Basic instantiation of a geometric representation of a molecule. Only include
    methods pertaining to what's possible with a 2D representation. Any methods that
    would result in a more defined representation (e.g. 3D, MD optimized, DFT optimized)
    should yield the appropriate subclass.
    '''

    def __init__(self):
        self.contents = None

    # def load(self, path: str):
    #     """Load in the data file"""
    #     with open(path, 'r') as f:
    #         self.contents = f.readlines()
    #     self.path = path
    #     return self.contents

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

        # should return instance of Geometry3D
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


class Geometry3D(Geometry2D):
    '''
    Builds off of the 2D representation, with additional methods specifc to a
    representation with 3D coordinates. Any methods that would result in a more
    defined representation (e.g. MD optimized, DFT optimized) should yield
    the appropriate subclass.
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()

    def optimize(self, method='md', kwargs={}):
        '''Call appropriate optimization function and return result as the appropriate class'''
        if method.lower() == 'md':
            # don't actually return the below, just ensure the MD optimizer returns an
            # instance of that class
            return MDOptimizedGeometry()
        elif method.lower() == 'dft':
            # don't actually return the below, just ensure the DFT optimizer returns an
            # instance of that class
            return DFTOptimizedGeometry()
        else:
            raise ValueError('Optimization method "{}" not supported'.format(method))

        raise NotImplementedError


class MDOptimizedGeometry(Geometry3D):
    '''
    Builds off of the 3D representation, with additional methods specifc to a
    representation with MD optimized 3D coordinates. Any methods that would
    result in a more defined representation (e.g. DFT optimized) should yield
    the appropriate subclass.
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()


class DFTOptimizedGeometry(Geometry3D):
    '''
    Builds off of the 3D representation, with additional methods specifc to a
    representation with DFT optimized 3D coordinates. 
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()
