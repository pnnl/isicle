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
    if not check_path or (check_path and geom1.path == geom2.path):
        if not check_contents or (check_contents and geom1.contents == geom2.contents):
            if not check_mol or (check_mol and MolToSmiles(geom1.mol) == MolToSmiles(geom2.mol)):
                return True
    return False


@pytest.fixture()
def geom():
    return isicle.geometry.Geometry()


class TestLoad:

    @pytest.mark.parametrize('path,expected',
                             [('tests/resources/geom_test.pkl', ['C=C'])])
    def test_load_pickle(self, path, expected):

        # Initialize correctly saved pickle
        geom = isicle.geometry.load_pickle(path)

        # Test for expected smiles
        assert geom.contents == expected

    @pytest.mark.parametrize('path,expected',
                             [('tests/resources/geom_test.smi', IOError),
                              ('tests/resources/geom_test_bad.pkl', TypeError)])
    def test_load_pickle_fail(self, path, expected):
        with pytest.raises(expected):
            isicle.geometry.load_pickle(path)

    @ pytest.mark.parametrize('path,expected,saved_pkl',
                              [('tests/resources/geom_test.smi', ['C=C'], 'tests/resources/geom_test.pkl')])
    def test_load_smiles(self, path, expected, saved_pkl):

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
        assert compare(geom1, geom3)

    @ pytest.mark.parametrize('path,expected,saved_pkl',
                              [('tests/resources/geom_test.inchi', ['InChI=1S/C2H4/c1-2/h1-2H2'], 'tests/resources/geom_test.pkl')])
    def test_load_inchi(self, path, expected, saved_pkl):

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

    @ pytest.mark.parametrize('path,expected,saved_pkl',
                              [('tests/resources/geom_test.smarts', ['[#6]=[#6]'], 'tests/resources/geom_test.pkl')])
    def test_load_smarts(self, path, expected, saved_pkl):

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

    @ pytest.mark.parametrize('path,saved_pkl',
                              [('tests/resources/geom_test.mol', 'tests/resources/geom_test.pkl')])
    def test_load_mol(self, path, saved_pkl):

        # Initialize using direct call
        geom1 = isicle.geometry.load_mol(path)

        # Initialize using indirect call
        geom2 = isicle.geometry.load(path)

        # Test both routes ended in same place
        assert compare(geom1, geom2)

        # Note: not comparing to saved molecule because Mols returned from
        # load_mol do not have Hs explicitly added.

    @ pytest.mark.parametrize('path,saved_pkl',
                              [('tests/resources/geom_test.mol2', 'tests/resources/geom_test.pkl')])
    def test_load_mol2(self, path, saved_pkl):

        # Initialize using direct call
        geom1 = isicle.geometry.load_mol2(path)

        # Initialize using indirect call
        geom2 = isicle.geometry.load(path)

        # Test both routes ended in same place
        assert compare(geom1, geom2)

        # Note: not comparing to saved molecule because Mols returned from
        # load_mol2 do not have Hs explicitly added.

    @ pytest.mark.parametrize('path,saved_pkl',
                              [('tests/resources/geom_test.pdb', 'tests/resources/geom_test.pkl')])
    def test_load_pdb(self, path, saved_pkl):

        # Initialize using direct call
        geom1 = isicle.geometry.load_pdb(path)

        # Initialize using indirect call
        geom2 = isicle.geometry.load(path)

        # Test both routes ended in same place
        assert compare(geom1, geom2)

        # Note: not comparing to saved molecule because Mols returned from
        # load_mol2 do not have Hs explicitly added.

# class TestGeometry:
#
#     def test_init(self, geom):
#         assert isinstance(geom, isicle.geometry.Geometry)
#
#     # @pytest.mark.parametrize('path,expected',
#     #                          [('resources/geom_test.smi', 'resources/geom_test.mol'),
#     #                           ('resources/geom_test.inchi', 'resources/geom_test.mol'),
#     #                           ('resources/geom_test.xyz', 'resources/geom_test.mol')])
#
#     #TODO: implement test
#     def test_get_mol(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test__handle_inplace(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_desalt(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_neutralize(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_tautomerize(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_optimize(self, geom):
#         raise NotImplementedError
#
#     # TODO: additional tests for charged compounds
#     @pytest.mark.parametrize('path, expected',
#                               ['resources/geom_test_3D.mol', 0])
#     def test_total_partial_charge(self, geom, path, expected):
#         geom.load(localfile(path))
#         result = geom.total_partial_charge()
#
#         # test attribute
#         assert geom.result == expected
#
#         # test return
#         assert result == expected
#
#
#     @pytest.mark.parametrize('path, expected',
#                               ['resources/geom_test_3D.mol', 6])
#     def test_natoms(self, geom, path, expected):
#         geom.load(localfile(path))
#         result = geom.natoms()
#
#         # test attribute
#         assert geom.result == expected
#
#         # test return
#         assert result == expected
#
#     #TODO: implement test
#     def test_copy(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_to_smiles(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_to_inchi(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_to_smarts(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_to_xyz(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_to_pdb(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_save_smiles(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_save_inchi(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_save_smarts(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_save_xyz(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_save_pdb(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_save_mol(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_save_pickle(self, geom):
#         raise NotImplementedError
#
#     #TODO: implement test
#     def test_save(self, geom):
#         raise NotImplementedError

# TODO: implement test class for MDOptimizedGeometry
# TODO: implement test class for DFTOptimizedGeometry
