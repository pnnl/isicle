import pytest
import isicle
from isicle.adducts import _ionize_method_selector, load_ions, parse_ions, ionize, ExplicitIonizationWrapper, CRESTIonizationWrapper
from isicle.geometry import Geometry
import os


def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)


@pytest.fixture()
def explicit():
    return ExplicitIonizationWrapper()


@pytest.fixture()
def crest():
    return CRESTIonizationWrapper()


@pytest.mark.parametrize('method,expected',
                         [('explicit', ExplicitIonizationWrapper),
                          ('Explicit', ExplicitIonizationWrapper),
                          ('EXPLICIT', ExplicitIonizationWrapper),
                          ('CREST', CRESTIonizationWrapper),
                          ('crest', CRESTIonizationWrapper)])
def test__ionize_method_selector(method, expected):
    # Check correct class is yielded
    assert isinstance(_ionize_method_selector(method), expected)


@pytest.mark.parametrize('method',
                         [('xtb'),
                          ('expicit'),
                          ('openbabel')])
def test__ionize_method_selector_fail(method):
    # Check unsupported selections
    with pytest.raises(ValueError):
        _ionize_method_selector(method)


def test_ionize(explicit):
    # Load geometry externally
    geom = isicle.geometry.load(localfile('resources/geom_test.mol'))
    isicle.adducts.ionize(geom, ion_path=localfile(
        'resources/ion_list.txt'), ion_method='explicit')
