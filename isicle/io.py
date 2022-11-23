
import os
import pickle

import isicle
from rdkit import Chem


def _load_text(path: str):
    '''
    Load text from file.

    Parameters
    ----------
    path : str
        Path to text file.

    Returns
    -------
    list
        Lines from given text file.

    '''

    # Read file contents
    with open(path, 'r') as f:
        contents = f.readlines()

    # Strip each line
    return [x.strip() for x in contents]


def load_xyz(path):
    '''
    Load XYZ from file.

    Parameters
    ----------
    path : str
        Path to XYZ file.

    Return
    -------
    :obj:`~isicle.geometry.XYZGeometry`
        Molecule representation.

    '''

    # Initialize XYZGeometry instance
    geom = isicle.geometry.XYZGeometry()

    # Populate basename
    geom.basename = os.path.splitext(os.path.basename(path))[0]

    # Load xyz file contents
    geom.xyz = _load_text(path)

    return geom


def _check_mol(mol, string_struct):
    '''
    Check if mol failed to generate. If so, throw error.

    Parameters
    ----------
    mol : :obj:`~rdkit.Chem.rdchem.Mol'
        RDKit representation of molecule structure.
    string_struct : str
        Input used to initialize Mol object.

    '''

    if mol is None:
        raise ValueError(
            'Could not convert structure to mol: {}'.format(string_struct))


def _load_mol_from_file(path, func=None):
    '''
    Load RDKit mol representation from file (pdb, mol, mol2).

    Parameters
    ----------
    path : str
        Path to supported file.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    '''

    # Initialize geometry instance
    geom = isicle.geometry.Geometry()

    # Populate basename
    geom.basename = os.path.splitext(os.path.basename(path))[0]

    # Load mol representation
    mol = func(path, removeHs=False, strictParsing=False)
    _check_mol(mol, path)
    geom.mol = mol

    return geom


def load_mol(path):
    '''
    Load mol from file.

    Parameters
    ----------
    path : str
        Path to mol file.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    '''

    return _load_mol_from_file(path, func=Chem.MolFromMolFile)


def load_mol2(path: str):
    '''
    Load mol2 from file.

    Parameters
    ----------
    path : str
        Path to mol2 file.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    '''

    return _load_mol_from_file(path, func=Chem.MolFromMol2File)


def load_pdb(path):
    '''
    Load PDB from file.

    Parameters
    ----------
    path : str
        Path to PDB file.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    '''

    return _load_mol_from_file(path, func=Chem.MolFromPdbFile)


def _load_line_notation(path, func=None, force=False):
    '''
    Load line notation representation (InChI, SMILES, SMARTS) from file.

    Parameters
    ----------
    path : str
        Path to file
    force : bool
        Indicate whether to force load input, ignoring errors.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    '''

    # Initialize geometry instance
    geom = isicle.geometry.Geometry()

    # Populate basename
    geom.basename = os.path.splitext(os.path.basename(path))[0]

    # Load text
    text = _load_text(path)[0].strip()

    # Load without sanitization, downstream checks
    if force is True:
        mol = func(text, sanitize=False)
        _check_mol(mol, text)
        mol.UpdatePropertyCache(strict=False)

    # Safely load
    else:
        mol = func(text)
        _check_mol(mol, text)

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)
    _check_mol(mol, text)

    # Populate rdkit mol instance attribute
    geom.mol = mol

    return geom


def load_smiles(path, force=False):
    '''
    Load SMILES from file.

    Parameters
    ----------
    path : str
        Path to file.
    force : bool
        Indicate whether to force load input, ignoring errors.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    '''

    return _load_line_notation(path, func=Chem.MolFromSmiles, force=force)


def load_inchi(path, force=False):
    '''
    Load InChI from file.

    Parameters
    ----------
    path : str
        Path to file.
    force : bool
        Indicate whether to force load input, ignoring errors.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    '''

    return _load_line_notation(path, func=Chem.MolFromInchi, force=force)


def load_smarts(path, force=False):
    '''
    Load SMARTS from file.

    Parameters
    ----------
    path : str
        Path to file.
    force : bool
        Indicate whether to force load input, ignoring errors.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    '''

    return _load_line_notation(path, func=Chem.MolFromSmarts, force=force)


def load_pickle(path):
    '''
    Load pickled file.

    Parameters
    ----------
    path : str
        Path to pickle.

    Returns
    -------
    data
        Previously pickled object instance.

    '''

    # Load file
    with open(path, 'rb') as f:
        return pickle.loads(f)


def load(path, force=False):
    '''
    Reads in molecule information of the following supported file types:
    .smi, .inchi, .xyz, .mol, .mol2, .pkl, .pdb. Direct loaders can also
    be used, see load_* functions for more information.

    Parameters
    ----------
    path : str
        Path to file with molecule information.
    force : bool
        Indicate whether to force load input, ignoring errors.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry` or :obj:`~isicle.geometry.XYZGeometry`
        Molecule representation.

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
        return load_smiles(path, force=force)

    if extension == '.inchi':
        return load_inchi(path, force=force)

    if extension == '.smarts':
        return load_smarts(path)

    raise IOError('Extension {} not recognized.'.format(extension))
