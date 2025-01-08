import os
import pickle
from io import StringIO

import joblib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

import isicle


def _load_text(path: str):
    """
    Load text from file.

    Parameters
    ----------
    path : str
        Path to text file.

    Returns
    -------
    list
        Lines from given text file.

    """

    # Read file contents
    with open(path, "r") as f:
        contents = f.readlines()

    # Strip each line
    return [x.strip() for x in contents]


def load_xyz(path):
    """
    Load XYZ from file.

    Parameters
    ----------
    path : str
        Path to XYZ file.

    Return
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    """
    
    # # Check for charge
    # if charge is None:
    #     raise ValueError("Charge must be specified when loading XYZ files.")

    # Create mol object
    raw_mol = Chem.MolFromXYZFile(path)
    mol = Chem.Mol(raw_mol)
    # rdDetermineBonds.DetermineBonds(mol, charge=charge)

    # Basename
    basename = os.path.splitext(os.path.basename(path))[0]

    # Initialize Geometry instance
    geom = isicle.geometry.Geometry(mol=mol, basename=basename)

    return geom


def _check_mol(mol, string_struct):
    """
    Check if mol failed to generate. If so, throw error.

    Parameters
    ----------
    mol : :obj:`~rdkit.Chem.rdchem.Mol'
        RDKit representation of molecule structure.
    string_struct : str
        Input used to initialize Mol object.

    """

    if mol is None:
        raise ValueError("Could not convert structure to mol: {}".format(string_struct))


def _load_mol_from_file(path, func=None):
    """
    Load RDKit mol representation from file (pdb, mol, mol2).

    Parameters
    ----------
    path : str
        Path to supported file.
    func : onj
    basename : str

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    """

    # Load mol representation
    if func == Chem.MolFromMolFile:
        mol = func(path, removeHs=False, strictParsing=False)
    else:
        mol = func(path, removeHs=False)

    # Check result
    _check_mol(mol, path)

    # Populate basename
    basename = os.path.splitext(os.path.basename(path))[0]

    # Initialize geometry instance
    geom = isicle.geometry.Geometry(mol=mol, basename=basename)

    return geom


def load_mol(path):
    """
    Load mol from file.

    Parameters
    ----------
    path : str
        Path to mol file.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    """

    return _load_mol_from_file(path, func=Chem.MolFromMolFile)


def load_mol2(path: str):
    """
    Load mol2 from file.

    Parameters
    ----------
    path : str
        Path to mol2 file.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    """

    return _load_mol_from_file(path, func=Chem.MolFromMol2File)


def load_pdb(path):
    """
    Load PDB from file.

    Parameters
    ----------
    path : str
        Path to PDB file.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    """

    return _load_mol_from_file(path, func=Chem.MolFromPDBFile)


def _load_line_notation(path, func=None, force=False, string=False):
    """
    Load line notation representation (InChI, SMILES) from file.

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

    """

    basename = None
    if string:
        # Load text
        text = path

    else:
        # Load text
        text = _load_text(path)[0].strip()

        # Extract basename
        basename = os.path.splitext(os.path.basename(path))[0]

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

    if basename is None:
        basename = Chem.MolToInchiKey(mol)

    # Initialize geometry instance
    geom = isicle.geometry.Geometry(mol=mol, basename=basename)

    return geom


def load_smiles(path, force=False):
    """
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

    """

    extension = os.path.splitext(path)[-1].lower()
    if "smi" in extension:
        return _load_line_notation(
            path, func=Chem.MolFromSmiles, force=force
        )
    else:
        return _load_line_notation(
            path, func=Chem.MolFromSmiles, force=force, string=True
        )


def load_inchi(path, force=False):
    """
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

    """

    if "inchi=" in path.lower():
        return _load_line_notation(
            path, func=Chem.MolFromInchi, force=force, string=True
        )
    else:
        return _load_line_notation(
            path, func=Chem.MolFromInchi, force=force
        )


def load_pickle(path):
    """
    Load pickled file.

    Parameters
    ----------
    path : str
        Path to pickle.

    Returns
    -------
    data
        Previously pickled object instance.

    """

    # Load file
    with open(path, "rb") as f:
        return pickle.load(f)


def load_joblib(path):
    """
    Load joblib file.

    Parameters
    ----------
    path : str
        Path to pickle.

    Returns
    -------
    data
        Previously pickled object instance.

    """

    # Load file
    with open(path, "rb") as f:
        return joblib.load(f)


def _check_mol_obj(mol_obj):
    """ """
    
    if isinstance(mol_obj, Chem.Mol):
        return
    else:
        raise IOError("Not a valid RDKit Mol object passed.")


def load_mol_obj(mol_obj):
    """
    Load RDKit mol object into geometry instance

    Parameters
    ----------
    mol_obj : mol
        RDKit mol object

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    """

    # Validate mol object
    _check_mol_obj(mol_obj)

    # Initialize geometry instance
    geom = isicle.geometry.Geometry(mol=mol_obj,
                                    basename=Chem.MolToInchiKey(mol_obj))

    return geom


def load(path, **kwargs):
    """
    Reads in molecule information of the following supported file types:
    .smi, .inchi, .xyz, .mol, .mol2, .pkl, .pdb. Direct loaders can also
    be used, see load_* functions for more information.

    Parameters
    ----------
    path : str
        Path to file with molecule information.
    kwargs
        Keyword arguments passed to format-specific loaders.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry` or :obj:`~isicle.geometry.XYZGeometry`
        Molecule representation.

    """

    if isinstance(path, str):
        path = path.strip()
        extension = os.path.splitext(path)[-1].lower()

        if extension == ".pkl":
            return load_pickle(path)

        if "mol2" in extension:
            return load_mol2(path)

        if "mol" in extension:
            return load_mol(path)

        if extension == ".joblib":
            return load_joblib(path)

        if extension == ".xyz":
            return load_xyz(path, **kwargs)

        if extension == ".pdb":
            return load_pdb(path)

        if extension == ".inchi" or "inchi=" in path.lower():
            return load_inchi(path, **kwargs)

        try:
            return load_smiles(path, **kwargs)
        except:
            raise IOError("Extension {} not recognized.".format(extension))

    else:
        try:
            return load_mol_obj(path)
        except:
            raise IOError("Not a valid RDKit mol object passed.")


def save_xyz(path, geom):
    """
    Save molecule geometry as XYZ file.

    Parameters
    ----------
    path : str
        Path to output file.
    geom : :obj:`~isicle.geometry.Geometry`
            or :obj:`~isicle.geometry.XYZGeometry`
        Molecule representation.

    """

    # Check instance type
    if not isinstance(geom, isicle.geometry.Geometry):
        raise TypeError(
            "Must be `isicle.geometry.Geometry` to save in XYZ format."
        )

    # Write to file
    with open(path, "w") as f:
        f.write(geom.to_xyzblock())


def save_joblib(path, data):
    """
    Save object as joblib file.

    Parameters
    ----------
    path : str
        Path to output file.
    data : object
        Aribtrary object instance.

    """

    with open(path, "wb") as f:
        joblib.dump(data, f)


def save_pickle(path, data):
    """
    Save object as pickle file.

    Parameters
    ----------
    path : str
        Path to output file.
    data : object
        Aribtrary object instance.

    """

    with open(path, "wb") as f:
        pickle.dump(data, f)


def save_mfj(path, geom):
    """
    Save molecule geometry as MFJ file. Must have energy and charge information.

    Parameters
    ----------
    path : str
        Path to output file.
    geom : :obj:`~isicle.geometry.Geometry`
            or :obj:`~isicle.geometry.XYZGeometry`
        Molecule representation.

    """

    # Check instance type
    if not isinstance(geom, isicle.geometry.Geometry):
        raise TypeError(
            "Must be `isicle.geometry.Geometry` to save in XYZ format."
        )

    # Check for charges in global properties
    if (geom.energy is None) or (geom.charge is None):
        raise KeyError("DFT energy calculation required. See `isicle.qm.dft`.")

    # Get XYZ coordinates
    xyz = pd.read_csv(
        StringIO(geom.to_xyzblock()),
        skiprows=2,
        header=None,
        sep="\s+",
        names=["Atom", "x", "y", "z"],
    )

    # Extract and append charges
    xyz["Charge"] = geom.charge

    # Load masses and merge
    masses = isicle.utils.atomic_masses()[["Symbol", "Mass"]]
    mfj = pd.merge(xyz, masses, left_on="Atom", right_on="Symbol")

    # Rename columns
    mfj = mfj[["x", "y", "z", "Mass", "Charge"]].astype(float)

    # Write to file
    with open(path, "w") as f:
        f.write(os.path.splitext(os.path.basename(path))[0] + "\n")
        f.write("1\n")
        f.write(str(len(mfj.index)) + "\n")
        f.write("ang\n")
        f.write("calc\n")
        f.write("1.000\n")

        for row in mfj.values:
            f.write("\t".join([str(x) for x in row]) + "\n")


def save_smiles(path, geom):
    """
    Save molecule geometry as SMILES file.

    Parameters
    ----------
    path : str
        Path to output file.
    geom : :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    """

    # Check instance type
    if not isinstance(geom, isicle.geometry.Geometry):
        raise TypeError("Must be `isicle.geometry.Geometry` to save in SMILES format.")

    # Write
    with open(path, "w") as f:
        f.write(geom.to_smiles())


def save_inchi(path, geom):
    """
    Save molecule geometry as InChI file.

    Parameters
    ----------
    path : str
        Path to output file.
    geom : :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    """

    # Check instance type
    if not isinstance(geom, isicle.geometry.Geometry):
        raise TypeError("Must be `isicle.geometry.Geometry` to save in InChI format.")

    # Write
    with open(path, "w") as f:
        f.write(geom.to_inchi())


def save_mol(path, geom):
    """
    Save molecule geometry as MOL file.

    Parameters
    ----------
    path : str
        Path to output file.
    geom : :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    """

    # Check instance type
    if not isinstance(geom, isicle.geometry.Geometry):
        raise TypeError("Must be `isicle.geometry.Geometry` to save in MOL format.")

    # Write
    Chem.MolToMolFile(geom.mol, path)


def save_pdb(path, geom):
    """
    Save molecule geometry as PDB file.

    Parameters
    ----------
    path : str
        Path to output file.
    geom : :obj:`~isicle.geometry.Geometry`
        Molecule representation.

    """

    # Check instance type
    if not isinstance(geom, isicle.geometry.Geometry):
        raise TypeError("Must be `isicle.geometry.Geometry` to save in PDB format.")

    # Write
    Chem.MolToPDBFile(geom.mol, path)


def save(path, data):
    """
    Save molecule, format detected by path extension.

    Parameters
    ----------
    path : str
        Path to save file. Supported extensions include .pkl, .mfj, .xyz, .mol,
        .pdb, .inchi, .smi.
    data : obj
        Object instance. Must be :obj:`~isicle.geometry.Geometry` for .xyz and .mfj.

    """

    # Determine format from extension
    extension = os.path.splitext(path)[-1].lower()

    # Extension checks
    if extension == ".pkl":
        return save_pickle(path, data)

    if extension == ".joblib":
        return save_joblib(path, data)

    if extension == ".mfj":
        return save_mfj(path, data)

    if "mol" in extension:
        return save_mol(path, data)

    if extension == ".xyz":
        return save_xyz(path, data)

    if extension == ".pdb":
        return save_pdb(path, data)

    if "smi" in extension:
        return save_smiles(path, data)

    if extension == ".inchi":
        return save_inchi(path, data)

    raise IOError("Extension {} not recognized.".format(extension))
