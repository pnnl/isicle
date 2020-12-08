import pytest
import isicle
import os
import pandas as pd
import pickle

def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)

def compare(nwcr1, nwcr2):
    '''
    Compares two NWChemResult objects and returns if they are equivalent
    '''
    if \
    nwcr1.get_geometry().split('/')[-1] == nwcr2.get_geometry().split('/')[-1] and \
    nwcr1.get_energy() == nwcr2.get_energy() and \
    nwcr1.get_shielding() == nwcr2.get_shielding() and \
    nwcr1.get_spin() == nwcr2.get_spin() and \
    nwcr1.get_frequency() == nwcr2.get_frequency() and \
    nwcr1.get_molden() == nwcr2.get_molden():
        return True
    return False

@pytest.fixture()
def mparser():
    return isicle.parse.MobcalParser()

@pytest.fixture()
def iparser():
    return isicle.parse.ImpactParser()

@pytest.fixture()
def nparser():
    return isicle.parse.NWChemParser()


class TestMobcalParser:

    def test_init(self, mparser):
        assert isinstance(mparser, isicle.parse.MobcalParser)

    @pytest.mark.parametrize('path,expected',
                             [('resources/mobcal_output.txt', 241),
                              ('resources/mobcal_incomplete.txt', 220)])
    def test_load(self, mparser, path, expected):
        # initialize
        contents = mparser.load(localfile(path))

        # test attribute
        assert len(mparser.contents) == expected

        # test return
        assert len(contents) == expected

    @pytest.mark.parametrize('path,expected',
                             [('resources/mobcal_output.txt', {'ccs': [153.1950], 'std': [1.166770]}),
                              ('resources/mobcal_incomplete.txt', None)])
    def test_parse(self, mparser, path, expected):
        # initialize
        mparser.load(localfile(path))
        result = mparser.parse()

        # test attribute
        assert mparser.result == expected

        # test return
        assert result == expected

    # currently only tests success case
    @pytest.mark.parametrize('path,sep,nrows',
                             [('resources/mobcal_output.txt', '\t', 1)])
    def test_save(self, mparser, path, sep, nrows):
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

class TestImpactParser:

    def test_init(self, iparser):
        assert isinstance(iparser, isicle.parse.ImpactParser)

    @pytest.mark.parametrize('path,expected',
                             [('resources/impact_output.txt', 2)])
    def test_load(self, iparser, path, expected):
        # initialize
        contents = iparser.load(localfile(path))

        # test attribute
        assert len(iparser.contents) == expected

        # test return
        assert len(contents) == expected

    @pytest.mark.parametrize('path,expected_ccs,expected_semrel',
                             [('resources/impact_output.txt', [128.1511], [0.0013])])
    def test_parse(self, iparser, path, expected_ccs, expected_semrel):
        # initialize
        iparser.load(localfile(path))
        result = iparser.parse()

        # test attribute
        assert iparser.result['CCS_TJM'] == expected_ccs
        assert iparser.result['SEM_rel'] == expected_semrel

        # test return
        assert result['CCS_TJM'] == expected_ccs
        assert result['SEM_rel'] == expected_semrel

    # currently only tests success case
    @pytest.mark.parametrize('path,sep,nrows',
                             [('resources/impact_output.txt', '\t', 1)])
    def test_save(self, iparser, path, sep, nrows):
        # initialize
        output = localfile('resources/impact_save.txt')
        iparser.load(localfile(path))
        iparser.parse()
        iparser.save(output, sep=sep)

        # file exists
        assert os.path.exists(output)

        # read back in
        df = pd.read_csv(output, sep=sep)

        # check length
        assert len(df.index) == nrows

        # clean up
        os.remove(output)


class TestNWChemParser:

    def test_init(self, nparser):
        assert isinstance(nparser, isicle.parse.NWChemParser)

    @pytest.mark.parametrize('path,expected',
                             [('resources/nwchem_output/1R3R_difenacoum_+H_001_s.out', 444186)])
    def test_load(self, nparser, path, expected):

        # initialize
        contents = nparser.load(localfile(path))

        # test attribute
        assert len(nparser.contents) == expected

        # test return
        assert len(contents) == expected

    @pytest.mark.parametrize('path,expected_filename',
                             [('resources/nwchem_output/1R3R_difenacoum_+H_001_s.out',
                             'resources/nwchem_output/1R3R_difenacoum_+H_001_s.pkl')])
    def test_parse(self, nparser, path, expected_filename):
        # initialize
        nparser.load(localfile(path))
        result = nparser.parse(to_parse=['geometry', 'energy', 'shielding',
                                         'spin', 'frequency', 'molden'])

        # Load previous NWChemResult class
        expected = isicle.parse.NWChemResult()
        expected.load(expected_filename)

        # check against previously saved file
        assert compare(result, expected)

    # currently only tests success case
    @pytest.mark.parametrize('path',
                             [('resources/nwchem_output/1R3R_difenacoum_+H_001_s.out')])
    def test_save(self, nparser, path, temp_path='dummy.pkl'):

        # initialize
        nparser.load(localfile(path))
        result = nparser.parse(to_parse=['geometry', 'energy', 'shielding',
                                         'spin', 'frequency', 'molden'])
        nparser.save(temp_path)

        # file exists
        assert os.path.exists(temp_path)

        # read back in
        # Load previoud NWChemResult class
        saved_result = isicle.parse.NWChemResult()
        saved_result.load(temp_path)

        # check against newly saved file
        assert compare(result, saved_result)

        # clean up
        os.remove(temp_path)
