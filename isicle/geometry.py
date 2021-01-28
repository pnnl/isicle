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
# TODO: add docstrings


def _load_pickle(path):
    """Short summary.

    Parameters
    ----------
    path : type
        Description of parameter `path`.

    Returns
    -------
    type
        Description of returned object.

    """
    with open(path, 'rb') as f:
        mol = pickle.load(f)

    # Check for valid Geometry class type
    if mol.__class__.__name__ in ['Geometry', 'MDOptimizedGeometry',
                                  'DFTOptimizedGeometry']:
        return mol

    # This is not a Geometry* instance
    raise IOError('Unsupported geometry format: {}.'.format(mol.__class__.name))


def _load_text(path: str):
    '''Load in the data file'''
    with open(path, 'r') as f:
        contents = f.readlines()
    return contents


def _load_generic_geom(path):
    geom = Geometry()
    geom.path = path
    geom.contents = _load_text(path)
    geom.filetype = os.path.splitext(path)[-1].lower()
    return geom


def _load_xyz(path):
    geom = _load_generic_geom(path)
    xyz = next(pybel.readfile('xyz', path))
    name = ((path).split('.')[0]).split('/')[-1] + '.mol'
    geom.mol = xyz.write('mol', name, erwrite=True)
    return geom


def _load_mol(path):
    geom = _load_generic_geom(path)
    geom.mol = Chem.MolFromMolFile(path)
    return geom


# TODO: full implementation of pdb loader
def _load_pdb(path):
    geom = _load_generic_geom(path)
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

    if extension in ['mol', 'mol2']:
        return _load_mol(path)

    if extension == 'xyz':
        return _load_xyz(path)

    if extension == 'pdb':
        return _load_pdb(path)

    # No spatial info, return String instance
    if extension in ['smi', 'inchi', 'smiles']:
        return _load_generic_geom(path)

    raise IOError('Extension {} not recognized.'.format(extension))


# TODO: update documentation
class Geometry(GeometryInterface):
    '''
    Molecular representation as geometry with 3D coordinates.
    '''

    def __init__(self):
        self.path = None
        self.contents = None
        self.filetype = None
        self.mol = None

    def _gen_mol(self):
        """
        Uses loaded file with SMILES or InChI to generate a mol object and
        saves to mol attribute.

        Returns
        -------
        RDKit Mol Object
            Generated mol object
        """

        if self.filetype in ['smi', 'smiles', 'inchi']:
            self.mol = to_mol(self.contents, frm=self.filetype)
            return self.get_mol()  # Returns safe copy

        raise RuntimeError('mol object should have been generated upon load \
                           for the given filetype. Check original load() \
                           input.')

    # TODO: Return hard copy of mol object, not pointer
    def get_mol(self):
        return self.mol

    def desalt(self, salts=None, inplace=False):
        '''
        Converts SMILES to RDKit mol object.
        Desalts RDKit mol object using Chem.SaltRemover module.
        Returns desalted RDKit SMILES string.
        (not instituted) Salts to be removed with a config file.
        '''

        # Generate mol object if not already completed
        if self.mol is None:
            self._gen_mol()

        remover = SaltRemover.SaltRemover(defnFormat='smiles', defnData=salts)
        # defnData="[Cl,Br,Na]" *sample definition of salts to be removed*
        # add iterator for salts listed in config?
        # set potential salts to be removed in a config file

        mol, deleted = remover.StripMolWithDeleted(self.get_mol())
        # using StripMolWithDeleted instead of StripMol
        # add functionality to track removed salts
        # atomno = res.GetNumAtoms
        # if relevant to future use, returns atom count post desalting

        # Replace current smiles with desalted form
        if inplace:
            self.mol = mol

        return mol

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
            return [(to_mol(x, frm='smarts'), to_mol(y, frm='smi')) for x, y in patts]

        # Generate mol object if not already completed
        if self.mol is None:
            self._gen_mol()

        reactions = _InitializeNeutralisationReactions()

        mol = self.get_mol()
        replaced = False
        for i, (reactant, product) in enumerate(reactions):
            while mol.HasSubstructMatch(reactant):
                replaced = True
                rms = AllChem.ReplaceSubstructs(mol, reactant, product)
                mol = rms[0]

        # # TODO: is this still necessary?
        # if replaced:
        #     res = self.to_smiles(mol, frm='mol')
        # else:
        #     res = self.to_smiles(self.smiles, frm='smi')

        # Replace current smiles with neutral form
        if inplace:
            self.mol = mol

        return mol

    # TODO: Refactor based on new class structure
    def tautomerize(self, alt=False, inplace=False):
        '''
        Converts SMILES to RDKIT mol object.
        Generate tautomers according to RDKit TautomerEnumerator() method.
        Default returns first tautomer generated.
        If alt=True, returns all generated tautomers
        '''

        # Generate mol object if not already completed
        if self.mol is None:
            self._gen_mol()

        # source: https://rdkit.blogspot.com/2020/01/trying-out-new-tautomer.html
        # Discuss noted double bond changes
        enumerator = rdMolStandardize.TautomerEnumerator()

        mol = self.get_mol()
        res = [self.to_smiles()]
        tauts = enumerator.Enumerate(mol)
        smis = [to_smiles(x) for x in tauts]
        s_smis = sorted((x, y)
                        for x, y in zip(smis, tauts) if x != self.smiles)
        res += [y for x, y in s_smis]

        # Replace current smiles with major tautomer
        if inplace:
            self.smiles = res[0]  # TODO: change to update self.mol instead

        if alt is True:
            return res
        return res[0]

    def optimize(self, method='mtd', kwargs={}):

        # Generate mol object if not already completed
        if self.mol is None:
            self._gen_mol()

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

    def to_smiles(self):
        if self.mol is None:
            return to_smiles(self.contents, frm=self.filetype)
        return to_smiles(self.mol, frm='mol')

    def to_inchi(self):
        '''Return InChI string (in memory, no files).'''
        raise NotImplementedError

    def to_smarts(self):
        '''Return SMARTS string (in memory, no files).'''
        raise NotImplementedError

    def save_pickle(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        return

    def save_mol(self, path):
        '''Save RDKit Mol instance'''

        # Skip if no mol available
        if self.mol is None:
            self._gen_mol()

        # Save file
        with open(path, 'w+') as f:
            f.write(self.mol)
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


def smarts_to_mol(smarts):
    return Chem.MolFromSmarts(smarts)


def smi_to_mol(smi):
    mol = Chem.MolFromSmiles(smi)
    molH = Chem.AddHs(mol)
    mol = Chem.MolToMolBlock(molH)
    return mol


def inchi_to_mol(inchi):
    mol = Chem.MolFromInchi(inchi)
    molH = Chem.AddHs(mol)
    mol = Chem.MolToMolBlock(molH)
    return mol


# TODO: fix, need to generate xyz from contents in memory
def xyz_to_mol(path):
    '''Convert XYZ to Mol'''
    xyz = next(pybel.readfile('xyz', path))
    name = ((path).split('.')[0]).split('/')[-1] + '.mol'
    mol = xyz.write('mol', name, erwrite=True)
    return mol


def to_mol(cpd, frm='smi'):
    '''
    Convert from given datatype to RDKit Mol object.

    Can convert from smiles, smarts, inchi, or xyz.
    '''

    if frm.lower() in ['smiles', 'smi']:
        res = smi_to_mol(cpd)

    elif frm.lower() == 'smarts':
        res = smarts_to_mol(cpd)

    elif frm.lower() == 'inchi':
        res = inchi_to_mol(cpd)

    elif frm.lower() == 'xyz':
        res = xyz_to_mol(cpd)

    else:
        raise TypeError('Input format {} not supported.'.format(frm))

    if res is None:
        raise RuntimeError('Conversion to mol failed.')

    return res


def to_smiles(cpd, frm='mol'):
    '''
    Convert from given datatype to RDKit canonical SMILES.

    Can convert from SMILES, InChI or mol.
    '''

    if 'mol' in frm.lower():
        res = MolToSmiles(cpd)

    elif 'smi' in frm.lower():
        res = MolToSmiles(MolFromSmiles(cpd))

    elif 'inchi' in frm.lower():
        res = MolToSmiles(MolFromInchi(cpd))

    else:
        raise TypeError('Input format {} not supported.'.format(frm))

    if res is None:
        raise RuntimeError('Conversion to smiles failed.')

    return res


# TODO: still needed? What all should be supported?
def to_xyz(cpd, frm='mol'):
    '''
    Convert from given datatype to XYZ file.

    Can convert from smiles, inchi, or mol.
    '''

    if 'mol' in frm.lower():
        res = MolToXYZFile(cpd, filename)

    elif 'smi' in frm.lower():
        res = MolToXYZFile(smi_to_mol(cpd), filename)

    elif 'inchi' in frm.lower():
        res = MolToXYZFile(inchi_to_mol(cpd), filename)

    else:
        raise TypeError('Input format {} not supported.'.format(frm))

    if res is None:
        raise RuntimeError('Conversion to xyz failed.')

    return res
