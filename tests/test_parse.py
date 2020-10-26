import pytest
import isicle
import os
import pandas as pd


def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)


@pytest.fixture()
def mparser(path):
    return isicle.parse.MobcalParser()


class TestMobcalParser:
    @pytest.mark.parametrize('path',
                             ['resources/mobcal_output.txt',
                              'resources/mobcal_incomplete.txt'])
    def test_init(self, mparser, path):
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
