import isicle
import pytest
from rdkit import Chem
import os

from tests import localfile


@pytest.fixture
def geometry():
    geom = isicle.geometry.Geometry()
    geom.basename = 'test'
    mol = Chem.MolFromSmiles('CCCC')
    mol = Chem.AddHs(mol)
    Chem.AllChem.EmbedMolecule(mol)
    geom.mol = mol

    return geom


def test_load_xyz(geometry):
    # Path to file
    path = localfile('resources/test_load.xyz')

    # Save to file
    isicle.save(path, geometry)

    # Load
    geom = isicle.io.load_xyz(path)

    # Check instance is correct type
    assert isinstance(geom, isicle.geometry.XYZGeometry)

    # Check xyz attribute exists
    assert hasattr(geom, 'xyz')

    # Check xyz attribute is populated
    assert geom.xyz is not None

    # Clean up
    os.remove(path)


def test_load_mol(geometry):
    # Path to file
    path = localfile('resources/test_load.mol')

    # Save to file
    isicle.save(path, geometry)

    # Load
    geom = isicle.io.load_mol(path)

    # Check instance is correct type
    assert isinstance(geom, isicle.geometry.Geometry)

    # Check xyz attribute exists
    assert hasattr(geom, 'mol')

    # Check xyz attribute is populated
    assert geom.mol is not None

    # Clean up
    os.remove(path)


def test_load_mol2():
    # Path to file
    path = localfile('resources/geom_test.mol2')

    # Load
    geom = isicle.load(path)

    # Check instance is correct type
    assert isinstance(geom, isicle.geometry.Geometry)

    # Check xyz attribute exists
    assert hasattr(geom, 'mol')

    # Check xyz attribute is populated
    assert geom.mol is not None


def test_load_pdb(geometry):
    # Path to file
    path = localfile('resources/test_load.pdb')

    # Save to file
    isicle.save(path, geometry)

    # Load
    geom = isicle.io.load_pdb(path)

    # Check instance is correct type
    assert isinstance(geom, isicle.geometry.Geometry)

    # Check xyz attribute exists
    assert hasattr(geom, 'mol')

    # Check xyz attribute is populated
    assert geom.mol is not None

    # Clean up
    os.remove(path)


def test_load_smiles(geometry):
    # Path to file
    path = localfile('resources/test_load.smi')

    # Save to file
    isicle.save(path, geometry)

    # Load
    geom = isicle.io.load_smiles(path)

    # Check instance is correct type
    assert isinstance(geom, isicle.geometry.Geometry)

    # Check xyz attribute exists
    assert hasattr(geom, 'mol')

    # Check xyz attribute is populated
    assert geom.mol is not None

    # Clean up
    os.remove(path)


def test_load_inchi(geometry):
    # Path to file
    path = localfile('resources/test_load.inchi')

    # Save to file
    isicle.save(path, geometry)

    # Load
    geom = isicle.io.load_inchi(path)

    # Check instance is correct type
    assert isinstance(geom, isicle.geometry.Geometry)

    # Check xyz attribute exists
    assert hasattr(geom, 'mol')

    # Check xyz attribute is populated
    assert geom.mol is not None

    # Clean up
    os.remove(path)


def test_load_pickle(geometry):
    # Path to file
    path = localfile('resources/test_load.pkl')

    # Save to file
    isicle.save(path, geometry)

    # Load
    geom = isicle.io.load_pickle(path)

    # Check instance is correct type
    assert isinstance(geom, isicle.geometry.Geometry)

    # Check xyz attribute exists
    assert hasattr(geom, 'mol')

    # Check xyz attribute is populated
    assert geom.mol is not None

    # Clean up
    os.remove(path)


@pytest.mark.parametrize('path,instance',
                         [(localfile('resources/test_load.xyz'), isicle.geometry.XYZGeometry),
                          (localfile('resources/test_load.pkl'), isicle.geometry.Geometry),
                          # (localfile('resources/test_load.mfj'), isicle.geometry.XYZGeometry),
                          (localfile('resources/test_load.smi'), isicle.geometry.Geometry),
                          (localfile('resources/test_load.inchi'), isicle.geometry.Geometry),
                          (localfile('resources/test_load.mol'), isicle.geometry.Geometry),
                          # (localfile('resources/test_load.mol2'), isicle.geometry.Geometry),
                          (localfile('resources/test_load.pdb'), isicle.geometry.Geometry)
                         ])
def test_load(geometry, path, instance):
    # Save
    isicle.save(path, geometry)

    # Load
    geom = isicle.load(path)

    # Check not none
    assert geom is not None

    # Check correct type
    assert isinstance(geom, instance)

    # Clean up
    os.remove(path)


def test_save_xyz(geometry):
    # Output path
    path = localfile('resources/test_save.xyz')

    # Save to xyz
    isicle.io.save_xyz(path, geometry)

    # Check path exists
    assert os.path.exists(path)

    # Check not empty
    assert os.path.getsize(path) > 0

    # Clean up
    os.remove(path)


def test_save_pickle(geometry):
    # Output path
    path = localfile('resources/test_save.pkl')

    # Save to xyz
    isicle.io.save_pickle(path, geometry)

    # Check path exists
    assert os.path.exists(path)

    # Check not empty
    assert os.path.getsize(path) > 0

    # Clean up
    os.remove(path)


# def test_save_mfj():
#     # Output path
#     path = localfile('resources/test_save.mfj')

#     # Load optimized geometry
#     geom = isicle.load(localfile('resources/nwchem_output/methane.pkl'))

#     # Save
#     isicle.io.save_mfj(path, geom)

#     # Check path exists
#     assert os.path.exists(path)

#     # Check not empty
#     assert os.path.getsize(path) > 0

#     # Clean up
#     os.remove(path)


def test_save_smiles(geometry):
    # Output path
    path = localfile('resources/test_save.smi')

    # Save to xyz
    isicle.io.save_smiles(path, geometry)

    # Check path exists
    assert os.path.exists(path)

    # Check not empty
    assert os.path.getsize(path) > 0

    # Clean up
    os.remove(path)


def test_save_inchi(geometry):
    # Output path
    path = localfile('resources/test_save.inchi')

    # Save to xyz
    isicle.io.save_inchi(path, geometry)

    # Check path exists
    assert os.path.exists(path)

    # Check not empty
    assert os.path.getsize(path) > 0

    # Clean up
    os.remove(path)


def test_save_mol(geometry):
    # Output path
    path = localfile('resources/test_save.mol')

    # Save to xyz
    isicle.io.save_mol(path, geometry)

    # Check path exists
    assert os.path.exists(path)

    # Check not empty
    assert os.path.getsize(path) > 0

    # Clean up
    os.remove(path)


def test_save_pdb(geometry):
    # Output path
    path = localfile('resources/test_save.pdb')

    # Save to xyz
    isicle.io.save_pdb(path, geometry)

    # Check path exists
    assert os.path.exists(path)

    # Check not empty
    assert os.path.getsize(path) > 0

    # Clean up
    os.remove(path)


@pytest.mark.parametrize('path',
                         [(localfile('resources/test_save.xyz')),
                          (localfile('resources/test_save.pkl')),
                          # (localfile('resources/test_save.mfj')),
                          (localfile('resources/test_save.smi')),
                          (localfile('resources/test_save.inchi')),
                          (localfile('resources/test_save.mol')),
                          (localfile('resources/test_save.pdb'))
                         ])
def test_save(geometry, path):
    # Save
    isicle.save(path, geometry)

    # Check path exists
    assert os.path.exists(path)

    # Check not empty
    assert os.path.getsize(path) > 0

    # Clean up
    os.remove(path)
