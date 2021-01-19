from isicle.interfaces import GeometryGenerationInterface
from rdkit import Chem
from rdkit.Chem import AllChem
import pybel
import pickle
import numpy as np

# TO DO, read from xyz or mol2 to yield class
# in utils.py, read_mol, Mol, pop_aom, push_atom functions

# TODO: catch file not found errors

def _load_pickle(path):

    with open(path, 'rb') as f:
        mol = pickle.load(f)

    # Check for valid Geometry class type
    if instance.__class__.__name__ in ['Geometry', 'Geometry2D', 'Geometry3D']:
        return mol

    # This is not a Geometry* instance
    print('Error: unknown geometry loaded')
    return None

def _load_string_file(self, path: str):
    """Load in the data file"""
    with open(path, 'r') as f:
        contents = f.readlines()
    return contents

def _load_string(path):
    contents = _load_string_file(path)

    if contents is not None:
        Geom = Geometry()
        Geom.contents = contents
        Geom.path = path
        return Geom

def load(path, pkl=False):
    """
    Reads in molecule information of the below supported data types

    Supported file types: .smi, .inchi, .xyz, .mol, .mol2
    """
    if pkl or path.endswith('.pkl'):
        return _load_pickle(path)

    # Assume file with string information
    return _load_string(path)

# TODO: catch conversion error
def smarts_to_mol(smi):
    return Chem.MolFromSmarts(smi)

# TODO: catch conversion error
def smi_to_mol(smi):
    mol = Chem.MolFromSmiles(smi)
    molH = Chem.AddHs(mol)
    mol = Chem.MolToMolBlock(molH)
    return mol

# TODO: catch conversion error
def inchi_to_mol(inchi):
    mol = Chem.MolFromInchi(inchi)
    molH = Chem.AddHs(mol)
    mol = Chem.MolToMolBlock(molH)
    return mol

# TODO: catch conversion error
# TODO: fix, need to generate xyz from contents in memory
def xyz_to_mol(path):
    """Convert XYZ to Mol"""
    xyz = next(pybel.readfile("xyz", FILE))
    name = ((path).split('.')[0]).split('/')[-1] + '.mol'
    mol = xyz.write('mol', name, erwrite=True)
    return mol

def to_mol(self, text, from='SMILES'):
    """
    Convert from given datatype to RDKit Mol object.

    Can convert from smiles, smarts, inchi, or xyz.
    """

    if from.lower() in ['smiles', 'smi'] :
        return smi_to_mol(text)

    if from.lower() == 'smarts':
        return smarts_to_mol(text)

    if from.lower() == 'inchi':
        return inchi_to_mol(text)

    if from.lower() == 'xyz':
        return xyz_to_mol(text)

    # Failure
    print('Conversion failed. Input "from" type not supported.')
    return None

class Geometry(GeometryInterface):

    def __init__(self):
        self.path = None
        self.contents = None
        self.mol = None

    def to_2D(self):
        Geom2 = Geometry2D()
        Geom2.path = self.path
        Geom2.contents = self.contents
        Geom2.gen2D()
        return Geom2

    def to_3D(self):
        Geom3 = Geometry3D()
        Geom3.path = self.path
        Geom3.contents = self.contents
        Geom3.gen3D()
        return Geom3

    # TODO: update
    def total_partial_charge(self):
        if self.mol is None:
            return None
        self.charge = np.array([a.partialcharge for a in self.atoms]).sum()
        return self.charge

    # TODO: update
    def natoms(self):
        if self.mol is None:
            return None
        self.number = len(self.atoms)
        return self.number

    def save_pickle(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        return

    def save_mol(self, filename):
        if self.mol is None:
            print('A Mol object has not been generated yet.')
            return

        """Save RDKit Mol instance"""
        f = open(filename, 'w+')
        f.write(self.mol)
        f.close()
        return

    return


class Geometry2D(Geometry):
    '''
    Basic instantiation of a geometric representation of a molecule. Only include
    methods pertaining to what's possible with a 2D representation. Any methods that
    would result in a more defined representation (e.g. 3D, MD optimized, DFT optimized)
    should yield the appropriate subclass.
    '''

    def __init__(self):
        super().__init__()

    def gen2D(self, path: str):
        """Convert SMILES/INCHI into 2D geometry using RDKit"""
        try:
            extension = self.path.split('.')[-1]
        except IndexError:
            print('Invalid path')
            return None

        if extension in ['smi', 'inchi', 'xyz']:
            self.mol = to_mol(self.contents, from=extension)
            return self.mol

        elif self.path.endswith('.mol') or self.path.endswith('.mol2'):
            self.mol = self.contents
            return self.mol

        # Failure
        print('Invalid file extension')
        return None

    def save_mol_2D(self):
        self._save_mol('geom_2D.mol')
        return


# TODO: update documentation
class Geometry3D(Geometry):
    '''
    Builds off of the 2D representation, with additional methods specifc to a
    representation with 3D coordinates. Any methods that would result in a more
    defined representation (e.g. MD optimized, DFT optimized) should yield
    the appropriate subclass.
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()

    def gen3D(self, path: str, mol2D=None):
        """Convert 2D geometry mol file into 3D geometry using RDKit"""
        if mol2D is None:
            Geom2 = self.to_2D()  # Spawn dummy Geometry2D instance
            mol2D = Geom2.mol

        # Create 3D representation
        if mol2D is not None:
            self.mol = AllChem.MMFFOptimizeMolecule(mol2D)
            return self.mol

        # Failure
        return None

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

    def save_mol_3D(self):
        self.save_mol('geom_3D.mol')
        return

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
