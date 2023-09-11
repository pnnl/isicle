import os
import pickle
import joblib
from io import StringIO

import isicle
import pandas as pd
from rdkit import Chem
from openbabel import pybel


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


def load_xyz(path, basename=None):
    """
    Load XYZ from file.

    Parameters
    ----------
    path : str
        Path to XYZ file.
    basename : str
        Basename to save files as, if none default to inchikey

    Return
    -------
    :obj:`~isicle.geometry.XYZGeometry`
        Molecule representation.

    """

    # Initialize XYZGeometry instance
    geom = isicle.geometry.XYZGeometry()

    # Load xyz file contents
    geom.xyz = _load_text(path)

    # Create mol object
    mols = list(pybel.readfile("xyz", path))
    mo = mols[0]
    mo_block = mo.write("mol")
    geom.mol = Chem.MolFromMolBlock(mo_block, removeHs=False)

    # Populate basename
    if basename is None:
        try:
            geom.basename = os.path.splitext(os.path.basename(path))[0]
        except:
            geom.basename = Chem.MolToInchiKey(mol)
    else:
        geom.basename = basename

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


def _load_mol_from_file(path, func=None, basename=None):
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

    # Initialize geometry instance
    geom = isicle.geometry.Geometry()

    # Load mol representation
    if func == Chem.MolFromMolFile:
        mol = func(path, removeHs=False, strictParsing=False)
    else:
        mol = func(path, removeHs=False)

    # Check result
    _check_mol(mol, path)

    # Assign to instance attribute
    geom.mol = mol

    # Populate basename
    if basename is None:
        try:
            geom.basename = os.path.splitext(os.path.basename(path))[0]
        except:
            geom.basename = Chem.MolToInchiKey(mol)
    else:
        geom.basename = basename

    return geom


def load_nwchem(path, basename=None):
    """
    Load NWChem output from file.

    Parameters
    ----------
    path : str
        Path to nwchem output file

    Returns
    -------
    :obj:`~isicle.qm.NWChemWrapper`
        DFT representation

    """

    dft = isicle.qm.NWChemWrapper()

    parser = isicle.parse.NWChemParser()
    parser.load(path)
    result = parser.parse()

    dft.__dict__.update(result)
    dft.geom.add___dict__({k: v for k, v in result.items() if k != "geom"})
    dft.output = parser.load(path)

    if basename is None:
        try:
            basename = Chem.MolToInchiKey(dft.geom.mol)
        except:
            basename = os.path.splitext((os.path.basename(path)))[0]

    dft.basename = basename

    return dft


def load_mol(path, basename=None):
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

    return _load_mol_from_file(path, func=Chem.MolFromMolFile, basename=basename)


def load_mol2(path: str, basename=None):
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

    return _load_mol_from_file(path, func=Chem.MolFromMol2File, basename=basename)


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


def _load_line_notation(path, func=None, force=False, string=False, basename=None):
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

    # Initialize geometry instance
    geom = isicle.geometry.Geometry()

    if string:
        # Load text
        text = path
    else:
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

    # Populate basename
    if basename is None:
        geom.basename = Chem.MolToInchiKey(mol)
    else:
        geom.basename = basename

    return geom


def load_smiles(path, force=False, basename=None):
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
            path, func=Chem.MolFromSmiles, force=force, basename=basename
        )
    else:
        return _load_line_notation(
            path, func=Chem.MolFromSmiles, force=force, string=True, basename=basename
        )


def load_inchi(path, force=False, basename=None):
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
            path, func=Chem.MolFromInchi, force=force, string=True, basename=basename
        )
    else:
        return _load_line_notation(
            path, func=Chem.MolFromInchi, force=force, basename=basename
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
    try:
        Chem.MolToMolBlock(mol_obj)
    except ValueError:
        print("Invalid RDKit mol object")
        raise


def load_mol_obj(mol_obj, basename=None):
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

    # Initialize geometry instance
    geom = isicle.geometry.Geometry()

    # Validate mol object
    _check_mol_obj(mol_obj)

    # Populate rdkit mol instance attribute
    geom.mol = mol_obj

    # Populate basename
    if basename is None:
        geom.basename = Chem.MolToInchiKey(mol)
    else:
        geom.basename = basename

    return geom


def load(path, force=False, basename=None):
    """
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

    """
    if (type(path)) == str:
        path = path.strip()
        extension = os.path.splitext(path)[-1].lower()

        if extension == ".pkl":
            return load_pickle(path)

        if "mol2" in extension:
            return load_mol2(path, basename=basename)

        if "mol" in extension:
            return load_mol(path, basename=basename)

        if extension == ".joblib":
            return load_joblib(path)

        if extension == ".xyz":
            return load_xyz(path, basename=basename)

        if extension == ".pdb":
            return load_pdb(path, basename=basename)

        if extension == ".inchi" or "inchi=" in path.lower():
            return load_inchi(path, force=force, basename=basename)

        try:
            return load_smiles(path, force=force, basename=basename)
        except:
            raise IOError("Extension {} not recognized.".format(extension))

    else:
        try:
            return load_mol_obj(path, basename=basename)
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
    if not isinstance(geom, (isicle.geometry.Geometry, isicle.geometry.XYZGeometry)):
        raise TypeError(
            "Must be `isicle.geometry.Geometry` or \
                        `isicle.geometry.XYZGeometry` to save in XYZ format."
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
    if not isinstance(geom, (isicle.geometry.Geometry, isicle.geometry.XYZGeometry)):
        raise TypeError(
            "Must be `isicle.geometry.Geometry` or \
                        `isicle.geometry.XYZGeometry` to save in XYZ format."
        )

    # Check for charges in global properties
    if ("energy" not in geom.__dict__) or ("charge" not in geom.__dict__):
        raise KeyError("DFT energy calculation required. " "See isicle.qm.dft.")

    # Get XYZ coordinates
    xyz = pd.read_csv(
        StringIO(geom.to_xyzblock()),
        skiprows=2,
        header=None,
        delim_whitespace=True,
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


def save_sdf(path, geom):
    """
    Save molecule geometry(s) as SDF file.

    Parameters
    ----------
    path : str
        Path to output file.
    geom : :obj:`~isicle.geometry.Geometry` (or collection of)
        Molecule representation.
    """
    if isinstance(geom, isicle.geometry.Geometry):
        w = Chem.SDWriter(path)
        w.write(geom.mol)
        w = None

    elif isinstance(geom, isicle.geometry.Geometry):
        w = Chem.SDWriter(path)
        for m in geom.geom:
            w.write(m.mol)
        w = None


def save(path, data):
    """
    Save molecule, format detected by path extension.

    Parameters
    ----------
    path : str
        Path to save file. Supported extensions include .pkl, .mfj, .xyz, .mol,
        .pdb, .inchi, .smi.
    data : obj
        Object instance. Must be :obj:`~isicle.geometry.Geometry` or
        :obj:`~isicle.geometry.XYZGeometry` for .xyz and .mfj.

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
