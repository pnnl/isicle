from isicle.interfaces import GeometryInterface
from isicle.interfaces import MolecularStringInterface
from rdkit.Chem import SaltRemover, AllChem, MolToSmiles, MolFromSmiles, MolFromInchi, MolFromMolFile, MolToXYZFile
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem
import pybel
import pickle
import numpy as np
import os

# TO DO, read from xyz or mol2 to yield class
# in utils.py, read_mol, Mol, pop_aom, push_atom functions
# TODO: catch file not found errors
# TODO: add docstrings
# TODO: add relevant functions as class methods


def _load_pickle(path):
    with open(path, 'rb') as f:
        mol = pickle.load(f)

    # Check for valid Geometry class type
    if mol.__class__.__name__ in ['Geometry', 'MDOptimizedGeometry', 'DFTOptimizedGeometry']:
        return mol

    # This is not a Geometry* instance
    raise IOError('Unsupported geometry format: {}.'.format(mol.__class__.name))


def _load_text(path: str):
    '''Load in the data file'''
    with open(path, 'r') as f:
        contents = f.readlines()
    return contents


def _load_molecular_string(path):
    ms = MolecularString()
    ms.path = path
    ms.contents = _load_text(path)
    return ms


def _load_generic_3D(path):
    geom = Geometry()
    geom.path = path
    geom.contents = _load_text(path)

    return geom


def _load_xyz(path):
    geom = _load_generic_3D(path)
    xyz = next(pybel.readfile('xyz', path))
    name = ((path).split('.')[0]).split('/')[-1] + '.mol'
    geom.mol = xyz.write('mol', name, erwrite=True)
    return geom


def _load_mol(path):
    geom = _load_generic_3D(path)
    geom.mol = Chem.MolFromMolFile(path)
    return geom


# TODO: full implementation of pdb loader
def _load_pdb(path):
    geom = _load_generic_3D(path)
    return geom


def load(path):
    '''
    Reads in molecule information of the following supported file types:
    .smi, .inchi, .xyz, .mol, .mol2
    '''
    path = path.strip()
    extension = os.path.splitext(path)[-1].lower()

    if extension == 'pkl':
        return _load_pickle(path)

    elif extension in ['mol', 'mol2']:
        return _load_mol(path)

    elif extension == 'xyz':
        return _load_xyz(path)

    elif extension == 'pdb':
        return _load_pdb(path)

    # No spatial info, return String instance
    elif extension in ['smi', 'inchi']:
        return _load_molecular_string(path)

    raise IOError('Extension {} not recognized.'.format(extension))


class MolecularString(MolecularStringInterface):
    '''
    Refactor this. We know we either have InChI or SMILES (maybe SMARTS?)
    as input. Make use of this knowledge! If InChI, `to_inchi` just returns
    original string, `to_smiles` converts and returns. You could even just
    precompute all, such that the interfaces to other methods (e.g. `to_3D`)
    are the same by consistently using e.g. SMILES representation.
    '''

    def __init__(self):
        self.contents = None  # Original text from file
        self.path = None  # Path to original file
        self.smiles = None  # SMILES being manipulated

    def to_2D(self):
        # TODO: shouldn't rely on Geometry methods
        # We know we're starting with a string (inchi/smiles) here
        raise NotImplementedError

        # geom = Geometry()
        # geom.path = self.path
        # geom.contents = self.contents
        # return geom.to_2D()

    def to_3D(self):
        # TODO: shouldn't rely on Geometry methods
        # We know we're starting with a string (inchi/smiles) here
        raise NotImplementedError

        # geom = Geometry()
        # geom.path = self.path
        # geom.contents = self.contents
        # return geom.to_3D()

    def desalt(self, salts=None, inplace=False):
        '''
        Converts SMILES to RDKit mol object.
        Desalts RDKit mol object using Chem.SaltRemover module.
        Returns desalted RDKit SMILES string.
        (not instituted) Salts to be removed with a config file.
        '''
        remover = SaltRemover.SaltRemover(defnFormat='smiles', defnData=salts)
        # defnData="[Cl,Br,Na]" *sample definition of salts to be removed*
        # add iterator for salts listed in config?
        # set potential salts to be removed in a config file

        mol, deleted = remover.StripMolWithDeleted(self.to_mol(self.smiles))
        # using StripMolWithDeleted instead of StripMol
        # add functionality to track removed salts
        # atomno = res.GetNumAtoms
        # if relevant to future use, returns atom count post desalting

        res = self.to_smiles(mol, frm='mol')

        # Replace current smiles with desalted form
        if inplace:
            self.smiles = res

        return res

    def neutralize(self, inplace=False):
        '''
        Converts SMILES to RDKit mol object.
        Neutralizes SMILES string.
        Returns neutralized SMILES string.
        '''

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
            return [(self.to_mol(x, frm='smarts'), self.to_mol(y, frm='smi')) for x, y in patts]

        reactions = _InitializeNeutralisationReactions()

        mol = self.to_mol(self.smiles, frm='smi')
        replaced = False
        for i, (reactant, product) in enumerate(reactions):
            while mol.HasSubstructMatch(reactant):
                replaced = True
                rms = AllChem.ReplaceSubstructs(mol, reactant, product)
                mol = rms[0]
        if replaced:
            res = self.to_smiles(mol, frm='mol')
        else:
            res = self.to_smiles(self.smiles, frm='smi')

        # Replace current smiles with neutral form
        if inplace:
            self.smiles = res

        return res

    def tautomerize(self, alt=False, inplace=False):
        '''
        Converts SMILES to RDKIT mol object.
        Generate tautomers according to RDKit TautomerEnumerator() method.
        Default returns first tautomer generated.
        If alt=True, returns all generated tautomers
        '''
        # source: https://rdkit.blogspot.com/2020/01/trying-out-new-tautomer.html
        # Discuss noted double bond changes
        enumerator = rdMolStandardize.TautomerEnumerator()

        mol = self.to_mol(self.smiles, frm='smi')
        res = [self.smiles]
        tauts = enumerator.Enumerate(mol)
        smis = [self.to_smiles(x) for x in tauts]
        s_smis = sorted((x, y)
                        for x, y in zip(smis, tauts) if x != self.smiles)
        res += [y for x, y in s_smis]

        # Replace current smiles with major tautomer
        if inplace:
            self.smiles = res[0]

        if alt is True:
            return res
        else:
            return res[0]

    def to_smiles(self):
        '''Return SMILES string (in memory, no files).'''
        raise NotImplementedError

    def to_inchi(self):
        '''Return InChI string (in memory, no files).'''
        raise NotImplementedError

    def to_smarts(self):
        '''Return SMARTS string (in memory, no files).'''
        raise NotImplementedError

    def save(self, path, format=None):
        '''Generic save function, provided format flag.'''
        raise NotImplementedError


# TODO: update documentation
class Geometry(GeometryInterface):
    '''
    Molecular representation as geometry with 3D coordinates.
    '''

    def __init__(self):
        self.path = None
        self.contents = None
        self.mol = None

    def to_2D(self):
        # TODO: we have an internal mol representation, just
        # return using rdkit functionality
        raise NotImplementedError

    def optimize(self, method='mtd', kwargs={}):
        '''Call appropriate optimization function and return result as the appropriate class'''
        raise NotImplementedError

        if method.lower() == 'mtd':
            # don't actually return the below, just ensure the MD optimizer returns an
            # instance of that class
            return MDOptimizedGeometry()
        elif method.lower() == 'dft':
            # don't actually return the below, just ensure the DFT optimizer returns an
            # instance of that class
            return DFTOptimizedGeometry()
        else:
            raise ValueError(
                'Optimization method "{}" not supported'.format(method))

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
        return Chem.Mol.GetNumAtoms(self.mol)  # TODO: check

    def copy(self):
        cp = self.__class__
        cp.path = self.path
        cp.contents = self.contents
        cp.mol = self.mol  # TODO: make hard copy
        return cp

    def save_pickle(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        return

    def save_mol(self, path):
        '''Save RDKit Mol instance'''

        # Skip if no mol available
        if self.mol is None:
            print('A Mol object has not been generated yet.')
            return

        # Save file
        f = open(path, 'w+')
        f.write(self.mol)
        f.close()
        return

    def save(self, path, format=None):
        '''Generic save function, provided format flag.'''
        raise NotImplementedError


class MDOptimizedGeometry(Geometry):
    '''
    Builds off of the 3D representation, with additional methods specifc to a
    representation with MD optimized 3D coordinates. Any methods that would
    result in a more defined representation (e.g. DFT optimized) should yield
    the appropriate subclass.
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()


class DFTOptimizedGeometry(Geometry):
    '''
    Builds off of the 3D representation, with additional methods specific to a
    representation with DFT optimized 3D coordinates.
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()


# TODO: Roll all of the following into the functionality of MolecularString.
# Any converters we want to implement down the road will use the class, not the
# other way around (i.e. class uses functions).

# TODO: catch conversion error
def smarts_to_mol(path):
    return Chem.MolFromSmarts(path)


# TODO: catch conversion error
def smi_to_mol(path):
    init = load(path)
    mol = Chem.MolFromSmiles(init)
    molH = Chem.AddHs(mol)
    mol = Chem.MolToMolBlock(molH)
    return mol


# TODO: catch conversion error,
def inchi_to_mol(path):
    init = load(path)
    mol = Chem.MolFromInchi(init)
    molH = Chem.AddHs(mol)
    mol = Chem.MolToMolBlock(molH)
    return mol


# TODO: catch conversion error
# TODO: fix, need to generate xyz from contents in memory
def xyz_to_mol(path):
    '''Convert XYZ to Mol'''
    xyz = next(pybel.readfile('xyz', path))
    name = ((path).split('.')[0]).split('/')[-1] + '.mol'
    mol = xyz.write('mol', name, erwrite=True)
    return mol


def to_mol(path, frm='SMILES'):
    '''
    Convert from given datatype to RDKit Mol object.

    Can convert from smiles, smarts, inchi, or xyz.

    Need to read in smiles string from file, can't use filename as input
    '''

    if frm.lower() in ['smiles', 'smi']:
        return smi_to_mol(path)

    if frm.lower() == 'smarts':
        return smarts_to_mol(path)

    if frm.lower() == 'inchi':
        return inchi_to_mol(path)

    if frm.lower() == 'xyz':
        return xyz_to_mol(path)

    # Failure
    print('Conversion failed. Input "from" type not supported.')
    return None


def to_smiles(path, frm='mol'):
    '''
    Convert from given datatype to RDKit canonical SMILES.

    Can convert from smiles, inchi or mol.
    '''

    if frm.lower() == 'mol':
        return MolToSmiles(path)

    if frm.lower() == 'smiles':
        return MolToSmiles(MolFromSmiles(path))

    if frm.lower() == 'inchi':
        return MolToSmiles(MolFromInchi(path))

    # Failure
    print('Conversion failed. Input "from" type not supported.')
    return None


def to_xyz(path, filename, frm='mol'):
    '''
    Convert from given datatype to XYZ file.

    Can convert from smiles, inchi, or mol.
    '''

    if frm.lower() == 'mol':
        return MolToXYZFile(MolFromMolFile(path, sanitize=False, removeHs=False), filename)

    if frm.lower() == 'smiles':
        return MolToXYZFile(smi_to_mol(path), filename)

    if frm.lower() == 'inchi':
        return MolToXYZFile(inchi_to_mol(path), filename)

    # Failure
    print('Conversion failed. Input "from" type not supported.')
    return None
