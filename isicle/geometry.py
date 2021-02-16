import numpy as np
import os
import copy
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem
import pybel
import pickle
from isicle.interfaces import GeometryInterface


def load_pickle(path: str):
    '''
    Load pickled file.

    Parameters
    ----------
    path : type
        Path to pickle.

    Returns
    -------
    Geometry, MDOptimizedGeometry, or DFTOptimizedGeometry
        Previously pickled *Geometry instance.

    '''

    # Load file
    with open(path, 'rb') as f:
        try:
            mol = pickle.load(f)
        except pickle.UnpicklingError:
            raise IOError('Could not read file as pickle: {}'.format(path))

    # Check for valid Geometry class type
    if mol.__class__.__name__ in ['Geometry', 'MDOptimizedGeometry',
                                  'DFTOptimizedGeometry']:
        return mol

    # Failure. This is not a *Geometry instance
    raise TypeError('Unsupported geometry format: {}.'.format(mol.__class__))


def _load_text(path: str):
    '''
    Grab all text from given file.

    Parameters
    ----------
    path : str
        Path to text file.

    Returns
    -------
    list
        List of lines from given text file.

    '''
    with open(path, 'r') as f:
        contents = f.readlines()
    return [x.strip() for x in contents]


def _load_generic_geom(path: str):
    '''
    Create new Geometry instance and populate file information.

    Parameters
    ----------
    path : str
        Path to geometry text file (.mol, .smi, etc.)

    Returns
    -------
    Geometry
        Basic Geometry class with only file information populated.

    '''
    geom = Geometry()
    geom.path = path
    geom.contents = _load_text(path)
    geom.filetype = os.path.splitext(path)[-1].lower().strip()
    return geom


def load_xyz(path: str):
    '''
    Load XYZ file and return as a Geometry instance.

    Parameters
    ----------
    path : str
        Path to XYZ file

    Returns
    -------
    Geometry
        Provided file and molecule information

    '''
    # geom = _load_generic_geom(path)
    # xyz = next(pybel.readfile('xyz', path))

    # geom.mol = xyz.write('mol', None, overwrite=True)
    # return geom

    # NOTE: currently cannot cast to RDKit Mol object
    raise NotImplementedError


def load_mol(path: str):
    '''
    Load mol file and return as a Geometry instance.

    Parameters
    ----------
    path : str
        Path to mol file

    Returns
    -------
    Geometry
        Provided file and molecule information

    '''
    geom = _load_generic_geom(path)
    geom.mol = Chem.MolFromMolFile(path)
    return geom


def load_mol2(path: str):
    '''
    Load mol2 file and return as a Geometry instance.

    Parameters
    ----------
    path : str
        Path to mol2 file

    Returns
    -------
    Geometry
        Provided file and molecule information

    '''
    geom = _load_generic_geom(path)
    geom.mol = Chem.MolFromMol2File(path)
    return geom


def load_pdb(path: str):
    '''
    Load PDB file and return as a Geometry instance.

    Parameters
    ----------
    path : str
        Path to PDB file

    Returns
    -------
    Geometry
        Provided file and molecule information

    '''
    geom = _load_generic_geom(path)
    geom.mol = Chem.MolFromPDBFile(path)
    return geom


def _load_2D(path, convert_fxn):
    '''
    Load string file and return as a Geometry instance.

    Parameters
    ----------
    path : str
        Path to SMILES file
    convert_fxn: RDKit function
        Function to use to convert from string to mol (e.g. MolFromSmiles)

    Returns
    -------
    Geometry
        Provided file and molecule information

    '''
    geom = _load_generic_geom(path)
    string_struct = geom.contents[0].strip()
    mol = convert_fxn(string_struct)

    if mol is None:
        raise ValueError('Could not convert structure to mol: {}'.format(string_struct))

    # Hs not explicit, must be added. Not done for MolFromSmarts.
    if convert_fxn is not Chem.MolFromSmarts:
        mol = Chem.AddHs(mol)

    geom.mol = mol
    return geom


def load_smiles(path: str):
    '''
    Load SMILES file and return as a Geometry instance.

    Parameters
    ----------
    path : str
        Path to SMILES file

    Returns
    -------
    Geometry
        Provided file and molecule information

    '''
    return _load_2D(path, Chem.MolFromSmiles)


def load_inchi(path: str):
    '''
    Load InChI file and return as a Geometry instance.

    Parameters
    ----------
    path : str
        Path to InChI file

    Returns
    -------
    Geometry
        Provided file and molecule information

    '''
    return _load_2D(path, Chem.MolFromInchi)


def load_smarts(path: str):
    '''
    Load SMARTS file and return as a Geometry instance.

    Parameters
    ----------
    path : str
        Path to SMARTS file

    Returns
    -------
    Geometry
        Provided file and molecule information

    '''
    return _load_2D(path, Chem.MolFromSmarts)


def load(path: str):
    '''
    Reads in molecule information of the following supported file types:
    .smi, .inchi, .xyz, .mol, .mol2, .pkl, .pdb. Direct loaders can also
    be used, see load_* functions for more information.

    Parameters
    ----------
    path : str
        Path to file with molecule information.

    Returns
    -------
    Geometry
        Molecule information.

    '''
    path = path.strip()
    extension = os.path.splitext(path)[-1].lower()

    if extension == '.pkl':
        return load_pickle(path)

    if 'mol2' in extension:
        return load_mol2(path)

    if 'mol' in extension:
        return load_mol(path)

    if extension == '.xyz':
        return load_xyz(path)

    if extension == '.pdb':
        return load_pdb(path)

    if 'smi' in extension:
        return load_smiles(path)

    if extension == '.inchi':
        return load_inchi(path)

    if extension == '.smarts':
        return load_smarts(path)

    raise IOError('Extension {} not recognized.'.format(extension))


class Geometry(GeometryInterface):
    '''
    Molecule information, including information on the file it was
    generated from. It is not recommended to manipulate or retrieve
    attributes of this class without using class functions.

    Attributes
    ----------
    path : str
        Path provided to generate original instance.
    contents : list(str)
        Contents of file used to create original instance.
    filetype : str
        File type used to create original instance.
    mol : RDKit Mol object
        Current structure, potentially updated from its original
        form using functions in this class.

    '''

    def __init__(self, path=None, contents=None, filetype=None, mol=None):
        self.path = path
        self.contents = contents
        self.filetype = filetype
        self.mol = mol

    def get_mol(self, hard_copy=True):
        '''
        Returns RDKit Mol object for this Geometry.

        Parameters
        ----------
        hard_copy : boolean
            Return a hard copy of the mol object. If false, returns pointer to
            this instance's mol object (not recommended). Default: True.

        Returns
        -------
        RDKit Mol object
            Current structure

        '''
        return self.mol.__copy__()

    def _handle_inplace(self, mol, inplace):
        '''
        Return updated Geometry object with given structure.

        Parameters
        ----------
        mol : RDKit Mol object
            Structure to use.
        inplace : boolean
            If true, update this instance with the new structure. Otherwise,
            create a new Geometry instance and populate it with the structure.

        Returns
        -------
        Geometry
            Updated structure instance.

        '''
        if inplace:
            self.mol = mol
            return self

        # Make a new object and populate its mol with the given mol
        return type(self)(self.path, self.contents,
                          self.filetype, mol)

    def desalt(self, salts=None, inplace=False):
        '''
        Desalts RDKit mol object using Chem.SaltRemover module.

        Parameters
        ----------
        salts : str (optional)
            Salt type to remove. Ex: 'Cl', 'Br', '[Na+]'. Default: None.
        inplace : boolean (optional)
            If true, update this instance with the new structure. Otherwise,
            create a new Geometry instance and populate it with the structure.
            Default: False.

        Returns
        -------
        Geometry
            Molecule with given salt(S) removed.

        '''

        # If no salts given, skip desalting
        if salts is None:
            return self._handle_inplace(self.get_mol(), inplace)

        remover = SaltRemover(defnFormat='smiles', defnData=salts)
        # defnData="[Cl,Br,Na]" *sample definition of salts to be removed*
        # add iterator for salts listed in config?
        # set potential salts to be removed in a config file

        mol, deleted = remover.StripMolWithDeleted(self.get_mol())
        # using StripMolWithDeleted instead of StripMol
        # add functionality to track removed salts
        # atomno = res.GetNumAtoms
        # if relevant to future use, returns atom count post desalting

        return self._handle_inplace(mol, inplace)

    def neutralize(self, inplace=False):
        '''
        Neutralizes RDKit mol object using neutralization reactions.

        Parameters
        ----------
        inplace : boolean (optional)
            If true, update this instance with the new structure. Otherwise,
            create a new Geometry instance and populate it with the structure.
            Default: False.

        Returns
        -------
        Geometry
            Neutralized form of the molecule.

        '''

        def _initialize_neutralisation_reactions():
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
            return [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y)) for x, y in patts]

        reactions = _initialize_neutralisation_reactions()

        mol = self.get_mol()
        replaced = False
        for i, (reactant, product) in enumerate(reactions):
            while mol.HasSubstructMatch(reactant):
                replaced = True
                rms = Chem.AllChem.ReplaceSubstructs(mol, reactant, product)
                mol = rms[0]

        # # TODO: is this still necessary?
        # if replaced:
        #     res = self.to_smiles(mol, frm='mol')
        # else:
        #     res = self.to_smiles(self.smiles, frm='smi')

        return self._handle_inplace(mol, inplace)

    # TODO: Refactor based on new class structure
    def tautomerize(self, return_all=False, inplace=False):
        '''
        Generate tautomers according to RDKit TautomerEnumerator() method.

        Parameters
        ----------
        return_all : boolean (optional)
            If true, return all tautomers generated. Otherwise, only return
            the most common. Default=False
        inplace : boolean (optional)
            If true, update this instance with the new structure. Otherwise,
            create a new Geometry instance and populate it with the structure.
            Default: False.

        Returns
        -------
        Geometry or list(Geometry)
            Tautomer(s) of starting structure.

        '''

        # source: https://rdkit.blogspot.com/2020/01/trying-out-new-tautomer.html
        # Discuss noted double bond changes
        enumerator = rdMolStandardize.TautomerEnumerator()

        mol = self.get_mol()
        res = [mol]
        tauts = enumerator.Enumerate(mol)
        smis = [Chem.MolToSmiles(x) for x in tauts]
        s_smis = sorted((x, y)
                        for x, y in zip(smis, tauts) if x != self.to_smiles())
        res += [y for x, y in s_smis]

        # Ensure res is a list of mol objects
        if return_all:
            new_geoms = []
            for r in res:
                geom = self.__copy__()
                geom.mol = r
                new_geoms.append(geom)
            return new_geoms

        return self._handle_inplace(res[0], inplace)

    def dft_optimize(self, program='NWChem', template=None, **kwargs):
        '''
        Optimize geometry, either XYZ or PDB, using stated functional and basis set.
        Additional inputs can be grid size, optimization criteria level,
        '''
        # Select program
        qmw = _program_selector(program)

        # Load geometry
        qmw.set_geometry(self)

        # Save geometry
        qmw.save_geometry(path, fmt=kwargs.pop('fmt'))

        # Configure
        if template is not None:
            qmw.configure_from_template(template)
        else:
            qmw.configure(**kwargs)

        # Save configuration file
        qmw.save_config()

        # Run QM simulation
        qmw.run()

        # Finish/clean up
        return qmw.finish()

    # TODO: update
    def total_partial_charge(self):
        '''Sum the partial charge across all atoms.'''
        mol = self.get_mol()
        Chem.AllChem.ComputeGasteigerCharges(mol)
        contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge')
                    for i in range(mol.GetNumAtoms())]
        return np.nansum(contribs)

    def natoms(self):
        '''Calculate total number of atoms.'''
        return Chem.Mol.GetNumAtoms(self.get_mol())

    def __copy__(self):
        '''Return hard copy of this class instance.'''
        return type(self)(self.path, self.contents,
                          self.filetype, self.get_mol())

    def to_smiles(self):
        '''Get SMILES for this structure.'''
        return Chem.MolToSmiles(self.get_mol())

    def to_inchi(self):
        '''Get InChI for this structure.'''
        return Chem.MolToInchi(self.get_mol())

    def to_smarts(self):
        '''Get SMARTS for this structure.'''
        return Chem.MolToSmarts(self.get_mol())

    def to_xyzblock(self):
        #     '''Get XYZ text for this structure.'''
        #     return Chem.MolToXYZBlock(self.mol)
        # NOTE: Depricated, returns nothing for C2H4
        raise NotImplementedError

    def to_pdbblock(self):
        '''Get PDB text for this structure'''
        return Chem.MolToPDBBlock(self.get_mol())

    def to_molblock(self):
        '''Get PDB text for this structure'''
        return Chem.MolToMolBlock(self.get_mol())

    def save_smiles(self, path: str):
        '''Save this structure's SMILES to file.'''
        with open(path, 'w') as f:
            f.write(self.to_smiles())
        return 'Success'

    def save_inchi(self, path: str):
        '''Save this structure's InChI to file.'''
        with open(path, 'w') as f:
            f.write(self.to_inchi())
        return 'Success'

    def save_smarts(self, path: str):
        '''Save this structure's SMARTS to file.'''
        with open(path, 'w') as f:
            f.write(self.to_smarts())
        return 'Success'

    def save_xyz(self, path: str):
        #     '''Save XYZ file for this structure.'''
        #     return Chem.MolToXYZFile(self.get_mol(), path)
        # NOTE: Depricated, creates blank files for C2H4
        raise NotImplementedError

    def save_mol(self, path):
        '''Save Mol file for this structure.'''
        return Chem.MolToMolFile(self.get_mol(), path)

    def save_pickle(self, path):
        '''Pickle this class instance.'''
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        return 'Success'

    def save_pdb(self, path: str):
        '''Save PDB file for this structure.'''
        return Chem.MolToPDBFile(self.mol, path)

    def save(self, path, fmt=None):
        '''
        Save molecule

        Parameters
        ----------
        path : str
            Path to save file.
        fmt : str (optional)
            Format to save this molecule in. If None, determined from given
            path's extension. If .pkl. pickles this full object.
            Default: None.
            Supported formats: .smi (SMILES), .inchi (InChI), .mol, .xyz,
            .pdb, .pkl.

        Returns
        -------
        str
            Status of save.

        '''

        if fmt is None:
            # Decide format based on path
            fmt = os.path.splitext(path)[-1]
        fmt = fmt.lower()

        if 'smi' in fmt:
            return self.save_smiles(path)

        if 'inchi' in fmt:
            return self.save_inchi(path)

        if 'smarts' in fmt:
            return self.save_smarts(path)

        if 'xyz' in fmt:
            return self.save_xyz(path)

        if 'mol' in fmt:  # diff for mol2?
            return self.save_mol(path)

        if 'pkl' in fmt:
            return self.save_pickle(path)

        if 'pdb' in fmt:
            return self.save_pdb(path)

        # TODO: enable Compute2DCoords, https://www.rdkit.org/docs/source/rdkit.Chem.rdDepictor.html

        raise TypeError('Input format {} not supported.'.format(fmt))


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
