from isicle.interfaces import GeometryInterface
from isicle.interfaces import MolecularStringInterface
from rdkit.Chem import SaltRemover, AllChem, MolToSmiles, MolFromSmiles, MolToPDBFile
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
    if mol.__class__.__name__ in ['Geometry', 'Geometry2D', 'Geometry3D']:
        return mol

    # This is not a Geometry* instance
    # TODO: raise error
    print('Error: unknown geometry loaded')
    return None


def _load_text_file(self, path: str):
    """Load in the data file"""
    with open(path, 'r') as f:
        contents = f.readlines()
    return contents


def load(path, pkl=False):
    """
    Reads in molecule information of the below supported data types

    Supported file types: .smi, .inchi, .xyz, .mol, .mol2
    """
    path = path.strip()
    extension = os.path.splitext(path)[-1]

    if pkl or extension == 'pkl':
        return _load_pickle(path)

    # Spacial coordinates provided, return Geometry3D instance
    if extension in ['mol', 'mol2', 'xyz']:
        geom3 = Geometry3D()
        geom3.path = path
        geom3.contents = _load_text_file(path)

        if extension == 'xyz':
            xyz = next(pybel.readfile('xyz', FILE))
            name = ((self.path).split('.')[0]).split('/')[-1] + '.mol'
            geom3.mol = xyz.write('mol', name, erwrite=True)
        else:
            geom3.mol = Chem.MolFromMolFile(path)

        # Check for failure
        if geom3.mol is None:
            return None

        return geom3

    # No spatial info, return String instance
    if extension in ['smi', 'inchi']:
        ms = MolecularString()
        ms.path = path
        ms.contents = _load_text_file(path)
        ms.canon_smi = to_smi(ms.contents, frm=extension)
        return ms

    # Assume file with string information
    return _load_text(path)


# TODO: catch conversion error
def smarts_to_mol(smarts):
    return Chem.MolFromSmarts(smarts)


# TODO: catch conversion error
def smi_to_mol(smi):
    init = load(smi)
    mol = Chem.MolFromSmiles(init)
    molH = Chem.AddHs(mol)
    mol = Chem.MolToMolBlock(molH)
    return mol


# TODO: catch conversion error,
def inchi_to_mol(inchi):
    init = load(inchi)
    mol = Chem.MolFromInchi(init)
    molH = Chem.AddHs(mol)
    mol = Chem.MolToMolBlock(molH)
    return mol


# TODO: catch conversion error
# TODO: fix, need to generate xyz from contents in memory
def xyz_to_mol(path):
    """Convert XYZ to Mol"""
    xyz = next(pybel.readfile('xyz', path))
    name = ((path).split('.')[0]).split('/')[-1] + '.mol'
    mol = xyz.write('mol', name, erwrite=True)
    return mol


def to_mol(input, frm='SMILES'):
    """
    Convert from given datatype to RDKit Mol object.

    Can convert from smiles, smarts, inchi, or xyz.

    Need to read in smiles string from file, can't use filename as input
    """

    if frm.lower() in ['smiles', 'smi']:
        return smi_to_mol(input)

    if frm.lower() == 'smarts':
        return smarts_to_mol(input)

    if frm.lower() == 'inchi':
        return inchi_to_mol(input)

    if frm.lower() == 'xyz':
        return xyz_to_mol(input)

    # Failure
    print('Conversion failed. Input "from" type not supported.')
    return None


def to_smiles(input, frm='mol'):
    """
    Convert from given datatype to RDKit canonical SMILES.

    Can convert from smiles, inchi or mol.
    """

    if frm.lower() == 'mol':
        return MolToSmiles(input)

    if frm.lower() == 'smiles':
        return MolToSmiles(MolFromSmiles(input))

    if frm.lower() == 'inchi':
        return MolToSmiles(MolFromInchi(input))

    # Failure
    print('Conversion failed. Input "from" type not supported.')
    return None


def to_xyz(input, filename, frm='mol'):
    """
    Convert from given datatype to XYZ file.

    Can convert from smiles, inchi, or mol.
    """

    if frm.lower() == 'mol':
        return MolToXYZFile(MolFromMolFile(input, sanitize=False, removeHs=False), filename)

    if frm.lower() == 'smiles':
        return MolToXYZFile(smi_to_mol(input), filename)

    if frm.lower() == 'inchi':
        return MolToXYZFile(inchi_to_mol(input), filename)

    # Failure
    print('Conversion failed. Input "from" type not supported.')
    return None


class MolecularString(MolecularStringInterface):
    def __init__(self):
        self.contents = None  # Original text from file
        self.path = None  # Path to original file
        self.smiles = None  # SMILES being manipulated

    def to_Geometry(self, gen2D=False, gen3D=False):
        geom = Geometry()
        geom.path = self.path
        geom.contents = self.contents

        # Returns Geometry2D instance
        if gen2D:
            return geom.to_2D()

        # Returns Geometry3D instance
        if gen3D:
            return geom.to_3D()

        return geom

    def desalt(self, salts=None, inplace=False):
        """
        Converts SMILES to RDKit mol object.
        Desalts RDKit mol object using Chem.SaltRemover module.
        Returns desalted RDKit SMILES string.
        (not instituted) Salts to be removed with a config file.
        """
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
        """
        Converts SMILES to RDKit mol object.
        Neutralizes SMILES string.
        Returns neutralized SMILES string.
        """

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
        """
        Converts SMILES to RDKIT mol object.
        Generate tautomers according to RDKit TautomerEnumerator() method.
        Default returns first tautomer generated.
        If alt=True, returns all generated tautomers
        """
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

    def to_smiles(self, input, frm='smi'):
        return to_smiles(input, frm=frm)

    def to_mol(self, input, frm='smi'):
        return to_mol(input, frm=frm)

    def save_pdb(self, fn):
        """Write contents to string file"""
        with open(path, 'w') as f:
            f.write(self.smiles)
            f.close()

    # What is fn?
    def save_string(self, fn):
        """Write contents to PDB file"""
        MolToPDBFile(self.smiles, fn)


class Geometry(GeometryInterface):
    def __init__(self):
        self.path = None
        self.contents = None
        self.mol = None

    def _to_2D(self):
        try:
            extension = self.path.split('.')[-1]
        except IndexError:
            print('Invalid path')
            return None

        if extension in ['smi', 'inchi', 'xyz']:
            self.mol = to_mol(self.contents, frm=extension)
            return self.mol

        elif self.path.endswith('.mol') or self.path.endswith('.mol2'):
            self.mol = self.contents  # TODO: create mol object
            return self.mol

        # Failure
        print('Invalid file extension')
        return None

    def to_2D(self):
        if self.__class__.__name__ == 'Geometry2D':
            return self.copy()

        Geom2 = Geometry2D()
        Geom2.path = self.path
        Geom2.contents = self.contents
        Geom2.mol = self._to_2D()
        return Geom2

    def _to_3D(self, mol2D=None):
        """Convert 2D geometry mol file into 3D geometry using RDKit"""
        if mol2D is None:
            Geom2 = self.to_2D()  # Spawn dummy Geometry2D instance
            mol2D = Geom2.mol

        # Create 3D representation
        if mol2D is not None:
            self.mol = EmbedMolecule(MolFromMolFile(
                mol2D, sanitize=False, removeHs=False), randomSeed=0xf00d)
            return self.mol

        # TODO: Account for 3D coord already available - go straight to mol

        # Failure
        return None

    def to_3D(self):
        if self.__class__.__name__ == 'Geometry3D':
            return self.copy()

        Geom3 = Geometry3D()
        Geom3.path = self.path
        Geom3.contents = self.contents
        Geom3.mol = self._to_3D()
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
        return Chem.Mol.GetNumAtoms(self.mol)  # TODO: check

    def save_pickle(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        return

    def save_mol(self, filename):
        """Save RDKit Mol instance"""

        # Skip if no mol available
        if self.mol is None:
            print('A Mol object has not been generated yet.')
            return

        # Save file
        f = open(filename, 'w+')
        f.write(self.mol)
        f.close()
        return

    def copy(self):
        cp = self.__class__
        cp.path = self.path
        cp.contents = self.contents
        cp.mol = self.mol  # TODO: make hard copy
        return cp


class Geometry2D(Geometry):
    '''
    Basic instantiation of a geometric representation of a molecule. Only include
    methods pertaining to what's possible with a 2D representation. Any methods that
    would result in a more defined representation (e.g. 3D, MD optimized, DFT optimized)
    should yield the appropriate subclass.
    '''

    def __init__(self):
        super().__init__()

    def save_mol_2D(self):
        self.save_mol('geom_2D.mol')
        return


# TODO: update documentation
class Geometry3D(Geometry):
    '''
    Builds off of the 2D representation, with additional methods specifc to a
    representation with 3D coordinates. Any methods that would result in a more
    defined representation (e.g. MTD optimized, DFT optimized) should yield
    the appropriate subclass.
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()

    def conformer(self, method='mtd', kwarfs={}):
        if method.lower() == 'mtd':
            return MTD()

    def optimize(self, method='mtd', kwargs={}):
        '''Call appropriate optimization function and return result as the appropriate class'''
        if method.lower() == 'mtd':
            # don't actually return the below, just ensure the MD optimizer returns an
            # instance of that class
            return MTD()
        elif method.lower() == 'dft':
            # don't actually return the below, just ensure the DFT optimizer returns an
            # instance of that class
            return DFT()
        else:
            raise ValueError(
                'Optimization method "{}" not supported'.format(method))

        raise NotImplementedError

    def save_mol_3D(self):
        self.save_mol('geom_3D.mol')
        return


class MTD(Geometry3D):
    '''
    Builds off of the 3D representation, with additional methods specifc to a
    representation with MD optimized 3D coordinates. Any methods that would
    result in a more defined representation (e.g. DFT optimized) should yield
    the appropriate subclass.
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()


class DFT(Geometry3D):
    '''
    Builds off of the 3D representation, with additional methods specific to a
    representation with DFT optimized 3D coordinates.
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()
