import pytest
import isicle
import numpy as np
import os
from isicle.geometry import Geometry, MDOptimizedGeometry, DFTOptimizedGeometry


def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)


@pytest.fixture()
def random_values():
    rng = np.random.default_rng(12345)
    values = 1000 * rng.random(size=30)
    return values


@pytest.fixture()
def random_energies():
    rng = np.random.default_rng(54321)
    energies = rng.random(size=30)
    return energies


@pytest.fixture()
def conformers():
    x = [isicle.geometry.load(localfile('resources/geom_test.mol'))
         for i in range(30)]
    return isicle.conformers.ConformationalEnsemble(x)


@pytest.mark.parametrize('func,expected',
                         [('boltzmann', isicle.conformers.boltzmann),
                          ('simple', isicle.conformers.simple_average),
                          ('lowest', isicle.conformers.lowest_energy),
                          ('threshold', isicle.conformers.threshold)])
def test__function_selector(func, expected):
    # Check correct class is yielded
    assert isicle.conformers._function_selector(func) == expected


@pytest.mark.parametrize('func',
                         [('boltzman'),
                          ('smple'),
                          ('lowest-energy'),
                          ('thresh')])
def test__function_selector_fail(func):
    # Check unsupported selections
    with pytest.raises(ValueError):
        isicle.conformers._function_selector(func)


@pytest.mark.parametrize('func,expected',
                         [(isicle.conformers.boltzmann, True),
                          (isicle.conformers.simple_average, False),
                          (isicle.conformers.lowest_energy, True),
                          (isicle.conformers.threshold, True)])
def test__energy_based(func, expected):
    assert isicle.conformers._energy_based(func) is expected


def test_reduce(random_values, random_energies):
    result = isicle.conformers.reduce(random_values,
                                      func='boltzmann',
                                      energy=random_energies,
                                      index=None)

    assert abs(result['mean'] - 670.505) < 1E-3
    assert abs(result['std'] - 26.566) < 1E-3
    assert result['n'] == 30


@pytest.mark.parametrize('index',
                         [(None)])
def test_boltzmann(random_values, random_energies, index):
    result = isicle.conformers.boltzmann(random_values,
                                         random_energies,
                                         index=index)

    assert abs(result['mean'] - 670.505) < 1E-3
    assert abs(result['std'] - 26.566) < 1E-3
    assert result['n'] == 30


@pytest.mark.parametrize('index',
                         [(None)])
def test_simple_average(random_values, index):
    result = isicle.conformers.simple_average(random_values,
                                              index=index)

    assert abs(result['mean'] - 451.660) < 1E-3
    assert abs(result['std'] - 276.717) < 1E-3
    assert result['n'] == 30


@pytest.mark.parametrize('index',
                         [(None)])
def test_lowest_energy(random_values, random_energies, index):
    result = isicle.conformers.lowest_energy(random_values,
                                             random_energies,
                                             index=index)

    assert abs(result['value'] - 667.237) < 1E-3
    assert abs(result['energy'] - 0.016) < 1E-3


@pytest.mark.parametrize('index',
                         [(None)])
def test_threshold(random_values, random_energies, index):
    result = isicle.conformers.threshold(random_values,
                                         random_energies,
                                         threshold=0.5,
                                         index=index)

    assert abs(result['mean'] - 414.716) < 1E-3
    assert abs(result['std'] - 289.999) < 1E-3
    assert result['n'] == 15


@pytest.mark.parametrize('objects',
                         [([Geometry(), Geometry(), Geometry()]),
                          ([MDOptimizedGeometry(), DFTOptimizedGeometry(), Geometry()])])
def test_build_conformational_ensemble(objects):
    isicle.conformers.build_conformational_ensemble(objects)


@pytest.mark.parametrize('objects',
                         [([Geometry(), 'abc', Geometry()]),
                          ([MDOptimizedGeometry(), DFTOptimizedGeometry(), list()])])
def test_build_conformational_ensemble_fail(objects):
    with pytest.raises(TypeError):
        isicle.conformers.build_conformational_ensemble(objects)


class TestConformationalEnsemble:

    @pytest.mark.parametrize('index',
                             [(False)])
    def test_reduce(self, conformers, random_values, random_energies, index):
        # Set values
        for c, v, e in zip(conformers, random_values, random_energies):
            c.dummy = v
            c.energy = e

        # Reduce attribute
        result = conformers.reduce('dummy', func='boltzmann', index=index)

        # Verify result
        assert abs(result['mean'] - 670.505) < 1E-3
        assert abs(result['std'] - 26.566) < 1E-3
        assert result['n'] == 30

    def test__apply_method(self, conformers):
        result = conformers._apply_method('get_natoms')
        assert all(x == 2 for x in result)

    def test__apply_function(self, conformers):
        result = conformers._apply_function(isicle.qm.dft,
                                            program='NWChem',
                                            fmt='xyz')

    def test_apply(self, conformers):
        result = conformers.apply(method='get_natoms')
        assert all(x == 2 for x in result)
