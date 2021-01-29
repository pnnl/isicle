import pytest
import isicle
from isicle import geometry
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)


@pytest.fixture()
def geom():
    return isicle.geometry.Geometry()


class TestGeometry:

    def test_init(self, geom):
        assert isinstance(geom, isicle.geometry.Geometry)

    # @pytest.mark.parametrize('path,expected',
    #                          [('resources/geom_test.smi', 'resources/geom_test.mol'),
    #                           ('resources/geom_test.inchi', 'resources/geom_test.mol'),
    #                           ('resources/geom_test.xyz', 'resources/geom_test.mol')])

    #TODO: implement test
    def test_get_mol(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test__handle_inplace(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_desalt(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_neutralize(self, geom):
        raise NotImplementedError
    
    #TODO: implement test
    def test_tautomerize(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_optimize(self, geom):
        raise NotImplementedError

    # TODO: additional tests for charged compounds
    @pytest.mark.parametrize('path, expected',
                              ['resources/geom_test_3D.mol', 0])
    def test_total_partial_charge(self, geom, path, expected):
        geom.load(localfile(path))
        result = geom.total_partial_charge()
  
        # test attribute
        assert geom.result == expected

        # test return
        assert result == expected


    @pytest.mark.parametrize('path, expected',
                              ['resources/geom_test_3D.mol', 6])
    def test_natoms(self, geom, path, expected):
        geom.load(localfile(path))
        result = geom.natoms()

        # test attribute
        assert geom.result == expected

        # test return
        assert result == expected

    #TODO: implement test
    def test_copy(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_to_smiles(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_to_inchi(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_to_smarts(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_to_xyz(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_to_pdb(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_save_smiles(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_save_inchi(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_save_smarts(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_save_xyz(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_save_pdb(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_save_mol(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_save_pickle(self, geom):
        raise NotImplementedError

    #TODO: implement test
    def test_save(self, geom):
        raise NotImplementedError

# TODO: implement test class for MDOptimizedGeometry
# TODO: implement test class for DFTOptimizedGeometry
