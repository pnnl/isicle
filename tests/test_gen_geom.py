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


class TestGeometryGeneration:

    def test_init(self, geom):
        assert isinstance(geom, isicle.geometry.Geometry)

    @pytest.mark.parametrize('path,expected',
                             [('resources/geom_test.smi', 'resources/geom_test.mol'),
                              ('resources/geom_test.inchi', 'resources/geom_test.mol'),
                              ('resources/geom_test.xyz', 'resources/geom_test.mol')])
    def _test_to_2D(self, geom, path, expected):
        # initialize
        contents = geom._to_2D(localfile(path))
        # test attribute
        assert geom.contents == expected

        # test return
        assert contents == expected

    @pytest.mark.parametrize('path, expected',
                             [('resources/geom_test.smi', 'resources/geom_test.mol'),
                              ('resources/geom_test.xyz', 'resources/geom_test.mol'),
                              ('resources/geom_test.inchi', 'resources/geom_test.mol')])
    def test_to_2D(self, gengeom, path, expected):
        # initialize
        result = geom.to_2D(path)

        # test attribute
        assert geom.result == expected

        # test return
        assert result == expected

    @pytest.mark.parametrize('path, expected',
                             [('resources/geom_test.mol', 'resources/geom_test_3D.mol')])
    def _test_to_3D(self, gengeom, path, expected):
        # initialize
        result = geom._to_3D(path)

        # test attribute
        assert geom.result == expected

        # test return
        assert result == expected

    @pytest.mark.parametrize('path, expected',
                             [('resources/geom_test.mol', 'resources/geom_test_3D.mol')])
    def test_to_3D(self, gengeom, path, expected):
        # initialize
        result = geom.to_3D(path)

        # test attribute
        assert geom.result == expected

        # test return
        assert result == expected

    @pytest.mark.parametrize('path, expected',
                              ['resources/geom_test_3D.mol', 0])
    def test_total_partial_charge(self, gengeom, path, expected):
        gengeom.load(localfile(path))
        result = gengeom.total_partial_charge()
  
        # test attribute
        assert gengeom.result == expected

        # test return
        assert result == expected


    @pytest.mark.parametrize('path, expected',
                              ['resources/geom_test_3D.mol', 6])
    def test_natoms(self, gengeom, path, expected):
        gengeom.load(localfile(path))
        result = gengeom.natoms()

        # test attribute
        assert gengeom.result == expected

        # test return
        assert result == expected


