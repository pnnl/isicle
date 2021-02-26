import pytest
import isicle
from isicle import geometry
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, MolToSmiles


def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)


def compare(geom1, geom2, check_path=True, check_contents=True, check_mol=True):
    '''
    Compares two NWChemResult objects and returns if they are equivalent
    '''
    if not check_path or geom1.path == geom2.path:
        if not check_contents or geom1.contents == geom2.contents:
            if not check_mol or MolToSmiles(geom1.mol) == MolToSmiles(geom2.mol):
                return True
    return False


@pytest.fixture()
def xgeom():
    return isicle.geometry.load_xyz(localfile('resources/geom_test.xyz'))


@pytest.fixture()
def geom():
    return isicle.geometry.load_smiles(localfile('resources/geom_test.smi'))


@pytest.fixture()
def geom_salt():
    return isicle.geometry.load_smiles(localfile('resources/geom_test_salt.smi'))


@pytest.fixture()
def geom_taut():
    return isicle.geometry.load_smiles(localfile('resources/geom_test_taut.smi'))


class TestLoad:

    @pytest.mark.parametrize('path,expected',
                             [('resources/geom_test.pkl', ['C=C'])])
    def test_load_pickle(self, path, expected):

        # Initialize correctly saved pickle
        geom = isicle.geometry.load_pickle(localfile(path))

        # Test for expected smiles
        assert geom.contents == expected

    @pytest.mark.parametrize('path,expected',
                             [('resources/geom_test.smi', IOError),
                              ('resources/geom_test_bad.pkl', TypeError)])
    def test_load_pickle_fail(self, path, expected):
        with pytest.raises(expected):
            isicle.geometry.load_pickle(localfile(path))

    @ pytest.mark.parametrize('path,expected,saved_pkl',
                              [('resources/geom_test.smi', ['C=C'], 'resources/geom_test.pkl')])
    def test_load_smiles(self, path, expected, saved_pkl):

        path = localfile(path)
        saved_pkl = localfile(saved_pkl)

        # Initialize using direct call
        geom1 = isicle.geometry.load_smiles(path)

        # Initialize using indirect call
        geom2 = isicle.geometry.load(path)

        # Load saved result
        geom3 = isicle.geometry.load_pickle(saved_pkl)

        # Test for expected smiles
        assert geom1.contents == expected

        # Test both routes ended in same place
        assert compare(geom1, geom2)

        # Test this matches saved result
        assert compare(geom1, geom3, check_path=False)

    @ pytest.mark.parametrize('path,expected,saved_pkl',
                              [('resources/geom_test.inchi', ['InChI=1S/C2H4/c1-2/h1-2H2'], 'resources/geom_test.pkl')])
    def test_load_inchi(self, path, expected, saved_pkl):

        path = localfile(path)
        saved_pkl = localfile(saved_pkl)

        # Initialize using direct call
        geom1 = isicle.geometry.load_inchi(path)

        # Initialize using indirect call
        geom2 = isicle.geometry.load(path)

        # Load saved result
        geom3 = isicle.geometry.load_pickle(saved_pkl)

        # Test for expected inchi
        assert geom1.contents == expected

        # Test both routes ended in same place
        assert compare(geom1, geom2)

        # Test this matches saved result
        assert compare(geom1, geom3, check_path=False, check_contents=False)

    @ pytest.mark.parametrize('path,expected',
                              [('resources/geom_test.smarts', ['[#6]=[#6]'])])
    def test_load_smarts(self, path, expected):

        path = localfile(path)

        # Initialize using direct call
        geom1 = isicle.geometry.load_smarts(path)

        # Initialize using indirect call
        geom2 = isicle.geometry.load(path)

        # Test for expected inchi
        assert geom1.contents == expected

        # Test both routes ended in same place
        assert compare(geom1, geom2)

        # Note: not comparing to saved molecule because Mols returned from
        # SMARTS are different from those returned from other Mol gen fxns.

    @ pytest.mark.parametrize('path',
                              [('resources/geom_test.mol')])
    def test_load_mol(self, path):

        path = localfile(path)

        # Initialize using direct call
        geom1 = isicle.geometry.load_mol(path)

        # Initialize using indirect call
        geom2 = isicle.geometry.load(path)

        # Test both routes ended in same place
        assert compare(geom1, geom2)

        # Note: not comparing to saved molecule because Mols returned from
        # load_mol do not have Hs explicitly added.

    @ pytest.mark.parametrize('path',
                              [('resources/geom_test.mol2')])
    def test_load_mol2(self, path):

        path = localfile(path)

        # Initialize using direct call
        geom1 = isicle.geometry.load_mol2(path)

        # Initialize using indirect call
        geom2 = isicle.geometry.load(path)

        # Test both routes ended in same place
        assert compare(geom1, geom2)

        # Note: not comparing to saved molecule because Mols returned from
        # load_mol2 do not have Hs explicitly added.

    @ pytest.mark.parametrize('path',
                              [('resources/geom_test.pdb')])
    def test_load_pdb(self, path):

        path = localfile(path)

        # Initialize using direct call
        geom1 = isicle.geometry.load_pdb(path)

        # Initialize using indirect call
        geom2 = isicle.geometry.load(path)

        # Test both routes ended in same place
        assert compare(geom1, geom2)

        # Note: not comparing to saved molecule because Mols returned from
        # load_mol2 do not have Hs explicitly added.


class TestXYZGeometry:

    def test_init(self, xgeom):
        assert isinstance(xgeom, isicle.geometry.XYZGeometry)

    # TODO: implement test
    def test_dft_optimize(self, xgeom):
        raise NotImplementedError

    # TODO: implement test
    def test_md_optimize(self, xgeom):
        raise NotImplementedError

    def test_get_global_properties(self, xgeom):
        xgeom.calculate_global_properties()
        assert xgeom.global_properties == xgeom.get_global_properties()

    def test_get_natoms(self, xgeom):
        assert xgeom.get_natoms() == 6

    def test_calculate_global_properties(self, xgeom):
        d = xgeom.calculate_global_properties()
        assert d['natoms'] == 6

    def test_get_atom_indices(self, xgeom):

        # Test default (C & H)
        assert xgeom.get_atom_indices() == list(range(6))

        # Test H only
        assert xgeom.get_atom_indices(atoms=['H']) == list(range(2, 6))

    def test_to_xyzblock(self, xgeom):
        assert xgeom.to_xyzblock().split('\n') == xgeom.contents

    def test__copy__(self, xgeom):

        # Record original number of atoms
        # This is automatically stored in the obj's dict.
        starting_natoms = xgeom.get_natoms()

        # Make copy
        xgeom_cp = xgeom.__copy__()

        # Test copy is correct
        assert compare(xgeom, xgeom_cp, check_mol=False)

        # Test objects are not linked
        xgeom.global_properties['natoms'] = -1
        assert xgeom.get_global_properties()['natoms'] != starting_natoms
        assert xgeom_cp.get_global_properties()['natoms'] == starting_natoms

    def _test_save(self, geom, temp_path, expected=None):

        # Check file exists
        assert os.path.exists(temp_path)

        # Check smiles in file
        if expected is not None:
            assert isicle.geometry._load_text(temp_path) == expected

        # Remove temp file
        os.remove(temp_path)
#

    def test_save_xyz(self, xgeom, temp_path='temp.xyz'):

        # Test direct call
        xgeom.save_xyz(temp_path)
        self._test_save(xgeom, temp_path, None)

        # Test indirect callable
        xgeom.save(temp_path)
        self._test_save(xgeom, temp_path, None)

    def test_save_pickle(self, xgeom, temp_path='temp.pkl'):

        # Test direct call
        xgeom.save_pickle(temp_path)
        self._test_save(xgeom, temp_path, None)

        # Test indirect callable
        xgeom.save(temp_path)
        self._test_save(xgeom, temp_path, None)


class TestGeometry:

    def test_init(self, geom):
        assert isinstance(geom, isicle.geometry.Geometry)

    def test_get_mol(self, geom):
        assert isinstance(geom.mol, Chem.rdchem.Mol)

    @ pytest.mark.parametrize('expected',
                              [('O=[51Cr](=O)([O-])[O-]')])
    def test_desalt(self, geom_salt, expected):
        geom_salt.desalt(salts='[Na+]', inplace=True)
        assert geom_salt.to_smiles() == expected

    @ pytest.mark.parametrize('expected',
                              [('O=[51Cr](=O)([O-])[O-]')])
    def test__handle_inplace(self, geom_salt, expected):

        # Create desalted object
        geom_desalt = geom_salt.desalt(salts='[Na+]')

        # Make sure original Geometry not affected
        assert geom_salt.to_smiles() == 'O=[51Cr](=O)([O-])[O-].[Na+].[Na+]'

        # Test changing original in place
        geom_salt.desalt(salts='[Na+]', inplace=True)
        assert geom_salt.to_smiles() == expected

    @ pytest.mark.parametrize('expected',
                              [('O=[51Cr](=O)(O)O')])
    def test_neutralize(self, geom_salt, expected):
        geom_salt.desalt(salts='[Na+]', inplace=True)
        geom_salt.neutralize(inplace=True)
        assert geom_salt.to_smiles() == expected

    # TODO: test for tautomers diff from starting structure
    @ pytest.mark.parametrize('expected',
                              [('[H]C([H])([H])C(=O)[O-]', '[H]C([H])([H])=C([O-])O', '[H][CH]([H])([H])C(=O)[O-]')])
    def test_tautomerize(self, geom_taut, expected):
        tauts = geom_taut.tautomerize(return_all=True)
        taut_smis = [x.to_smiles() for x in tauts]
        assert set(expected) == set(taut_smis)

    # TODO: implement test
    def test_dft_optimize(self, geom):
        raise NotImplementedError

    # TODO: implement test
    def test_md_optimize(self, geom):
        raise NotImplementedError

    def test_get_total_partial_charge(self, geom_salt):
        assert geom_salt.get_total_partial_charge() == 2

    def test_get_natoms(self, geom_salt):
        assert geom_salt.get_natoms() == 7

    def test_get_atom_indices(self, geom):

        # Test default (C & H)
        assert geom.get_atom_indices() == list(range(6))

        # Test H only
        assert geom.get_atom_indices(atoms=['H']) == list(range(2, 6))

    def test__copy__(self, geom_salt):

        # Record original smiles
        starting_smi = geom_salt.to_smiles()

        # Make copy
        geom_cp = geom_salt.__copy__()

        # Test copy is correct
        assert compare(geom_salt, geom_cp)

        # Test objects are not linked
        geom_salt.desalt(salts='[Na+]', inplace=True)
        assert geom_salt.to_smiles() != starting_smi
        assert geom_cp.to_smiles() == starting_smi

    def test_to_smiles(self, geom):
        assert geom.to_smiles() == '[H]C([H])=C([H])[H]'

    def test_to_inchi(self, geom):
        assert geom.to_inchi() == 'InChI=1S/C2H4/c1-2/h1-2H2'

    def test_to_smarts(self, geom):
        assert geom.to_smarts() == '[#6](=[#6](-[H])-[H])(-[H])-[H]'

    def _test_save(self, geom, temp_path, expected=None):

        # Check file exists
        assert os.path.exists(temp_path)

        # Check smiles in file
        if expected is not None:
            assert isicle.geometry._load_text(temp_path) == expected

        # Remove temp file
        os.remove(temp_path)

    def test_save_smiles(self, geom, temp_path='temp.smi'):

        # Test direct call
        geom.save_smiles(temp_path)
        self._test_save(geom, temp_path, [geom.to_smiles()])

        # Test indirect callable
        geom.save(temp_path)
        self._test_save(geom, temp_path, [geom.to_smiles()])

    def test_save_inchi(self, geom, temp_path='temp.inchi'):

        # Test direct call
        geom.save_inchi(temp_path)
        self._test_save(geom, temp_path, [geom.to_inchi()])

        # Test indirect callable
        geom.save(temp_path)
        self._test_save(geom, temp_path, [geom.to_inchi()])

    def test_save_smarts(self, geom, temp_path='temp.smarts'):

        # Test direct call
        geom.save_smarts(temp_path)
        self._test_save(geom, temp_path, [geom.to_smarts()])

        # Test indirect callable
        geom.save(temp_path)
        self._test_save(geom, temp_path, [geom.to_smarts()])

    def test_save_pickle(self, geom, temp_path='temp.pkl'):

        # Test direct call
        geom.save_pickle(temp_path)
        self._test_save(geom, temp_path, None)

        # Test indirect callable
        geom.save(temp_path)
        self._test_save(geom, temp_path, None)

# TODO: implement test class for MDOptimizedGeometry
# TODO: implement test class for DFTOptimizedGeometry
