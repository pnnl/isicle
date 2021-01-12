import pytest
import isicle
import os
import pandas as pd


def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)


@pytest.fixture()
def gengeom():
    return isicle.generate_geom.GeometryGeneration()


class TestGeometryGeneration:

    def test_init(self, gengeom):
        assert isinstance(gengeom, isicle.generate_geom.GeometryGeneration)

    @pytest.mark.parametrize('path,expected',
                             [('resources/geom_test.smi', 23),
                              ('resources/geom_test.inchi', 25),
                              ('resources/geom_test.xyz', 23)])
    def test_load(self, gengeom, path, expected):
        # initialize
        contents = gengeom.load(localfile(path))
        print(contents)
        # test attribute
        assert len(gengeom.contents) == expected

        # test return
        assert len(contents) == expected

    @pytest.mark.parametrize('path, expected',
                             [('resources/geom_test.smi', 'resources/geom_2D_smi.mol'),
                              ('resources/geom_test.xyz', 'resources/geom_2D.mol'),
                              ('resources/geom_test.inchi', 'resources/geom_2D.mol')])
    def test_convert2D(self, gengeom, path, expected):
        # initialize
        gengeom.load(localfile(path))
        result = gengeom.inputto2D()

        # test attribute
        assert gengeom.result == expected

        # test return
        assert result == expected

    @pytest.mark.parametrize('path, expected',
                             [('resources/geom_2D.mol', 'resources/geom_3D_xyz.mol'),
                              ('resources/geom_2D_smi.mol', 'resources/geom_3D_smi.mol'),
                              ('resources/geom_2D_inchi.mol', 'resources/geom_3D_inchi.mol')])
    def test_convert3D(self, gengeom, path, expected):
        # initialize
        gengeom.load(localfile(path))
        result = gengeom.convert3D()

        # test attribute
        assert gengeom.result == expected

        # test return
        assert result == expected

    # currently only tests success case
    @pytest.mark.parametrize('path,sep,nrows',
                             [('resources/geom_output.mol', '\t', 1)])
    def test_save(self, gengeom, path, sep, nrows):
        # initialize
        output = localfile('resources/mobcal_save.txt')
        mparser.load(localfile(path))
        mparser.parse()
        mparser.save(output, sep=sep)

        # file exists
        assert os.path.exists(output)

        # read back in
        df = pd.read_csv(output, sep=sep)

        # check length
        assert len(df.index) == nrows

        # clean up
        os.remove(output)
