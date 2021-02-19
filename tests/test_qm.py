import pytest
import isicle
from isicle.qm import _program_selector, NWChemWrapper, dft
from isicle.geometry import Geometry
import os
import shutil


def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)


@pytest.fixture()
def nwc():
    return NWChemWrapper()

@pytest.mark.parametrize('program,expected',
                         [('NWChem', NWChemWrapper),
                          ('nwchem', NWChemWrapper),
                          ('NWCHEM', NWChemWrapper)])
def test__program_selector(program, expected):
    # Check correct class is yielded
    assert isinstance(_program_selector(program), expected)


@pytest.mark.parametrize('program',
                         [('nw'),
                          ('gaussian'),
                          ('orca')])
def test__program_selector_fail(program):
    # Check unsupported selections
    with pytest.raises(ValueError):
        _program_selector(program)


def test_dft(nwc):
    # Load geometry externally
    geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

    # Run dft
    result = dft(geom, program='NWChem', fmt='xyz')


class TestNWChemWrapper:

    def test_init(self, nwc):
        # Check class initialization
        assert isinstance(nwc, NWChemWrapper)

        # Check temp directory created
        assert os.path.exists(nwc.temp_dir.name)

        # Clean up
        nwc.temp_dir.cleanup()

    @pytest.mark.parametrize('path,expected',
                             [('resources/geom_test.mol', 'geom_test')])
    def test_set_geometry(self, nwc, path, expected):
        # Load geometry externally
        geom = isicle.geometry.load(localfile(path))

        # Set geometry
        nwc.set_geometry(geom)

        # Check class instance
        assert isinstance(nwc.geom, Geometry)

        # Check basename attribute
        assert nwc.geom.basename == expected

        # Clean up
        nwc.temp_dir.cleanup()

    @pytest.mark.parametrize('fmt',
                             [('xyz'),
                              ('pdb')])
    def test_save_geometry(self, nwc, fmt):
        # Load geometry externally
        geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

        # Set geometry
        nwc.set_geometry(geom)

        # Save geometry
        nwc.save_geometry(fmt=fmt)

        # Check if file exists
        assert os.path.exists(os.path.join(nwc.temp_dir.name,
                                           '{}.{}'.format(nwc.geom.basename,
                                                          fmt.lower())))

        # Clean up
        nwc.temp_dir.cleanup()

    # TODO: parametrize
    # TODO: meaningful assertions (currently just runs without error)
    @pytest.mark.parametrize('tasks,ao_basis,cosmo',
                             [('optimize', 'spherical', True),
                              ('shielding', 'spherical', False),
                              ('spin', 'cartesian', True),
                              (['optimize', 'shielding'], ['spherical', 'spherical'], [False, True]),
                              (['shielding', 'spin'], 'spherical', True),
                              (['optimize', 'shielding', 'spin'], 'spherical', True)])
    def test_configure(self, nwc, tasks, ao_basis, cosmo):
        # Load geometry externally
        geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

        # Set geometry
        nwc.set_geometry(geom)

        # Save geometry
        nwc.save_geometry(fmt='pdb')

        # Configure
        config = nwc.configure(tasks=tasks, ao_basis=ao_basis, cosmo=cosmo)

        # Clean up
        nwc.temp_dir.cleanup()

    @pytest.mark.parametrize('basename,dirname',
                             [(None, None),
                              (None, 'override'),
                              ('override', None),
                              ('override', 'override')])
    def test_configure_from_template(self, nwc, basename, dirname):
        # Load geometry externally
        geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

        # Set geometry
        nwc.set_geometry(geom)

        # Save geometry
        nwc.save_geometry(fmt='pdb')

        # Configure from template
        config = nwc.configure_from_template(localfile('resources/nwchem_template.txt'),
                                             basename_override=basename,
                                             dirname_override=dirname,
                                             charge=0)

        # Clean up
        nwc.temp_dir.cleanup()

    def test_save_config(self, nwc):
        # Load geometry externally
        geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

        # Set geometry
        nwc.set_geometry(geom)

        # Save geometry
        nwc.save_geometry(fmt='pdb')

        # Configure
        nwc.configure()

        # Save config
        nwc.save_config()

        assert os.path.exists(os.path.join(nwc.temp_dir.name,
                                           nwc.geom.basename + '.nw'))

        # Clean up
        nwc.temp_dir.cleanup()

    # Not particularly testable
    def test_run(self, nwc):
        # Load geometry externally
        geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

        # Set geometry
        nwc.set_geometry(geom)

        # Save geometry
        nwc.save_geometry(fmt='pdb')

        # Configure
        nwc.configure()

        # Save config
        nwc.save_config()

        # For now, just check if NWChem is added to path and returns
        nwc.run()

        # Clean up
        nwc.temp_dir.cleanup()

    def test_finish(self, nwc):
        # Load geometry externally
        geom = isicle.geometry.load(localfile('resources/nwchem_output/1R3R_difenacoum_+H_001_s.xyz'))

        # Set geometry
        nwc.set_geometry(geom)

        # Copy example output to temp folder (i.e. assume nwchem was "run")
        shutil.copy2(localfile('resources/nwchem_output/1R3R_difenacoum_+H_001_s.out'),
                     nwc.temp_dir.name)

        # Finish
        nwc.finish(keep_files=False)
