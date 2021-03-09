import pytest
import isicle
from isicle.md import _program_selector, XTBWrapper, md
from isicle.geometry import Geometry
import os
import shutil


def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)


@pytest.fixture()
def xtb():
    return XTBWrapper()


@pytest.mark.parametrize('program,expected',
                         [('xtb', XTBWrapper),
                          ('XTB', XTBWrapper),
                          ('xTb', XTBWrapper)])
def test__program_selector(program, expected):
    # Check correct class is yielded
    assert isinstance(_program_selector(program), expected)


@pytest.mark.parametrize('program',
                         [('xt'),
                          ('amber')])
def test__program_selector_fail(program):
    # Check unsupported selections
    with pytest.raises(ValueError):
        _program_selector(program)


def test_md(xtb):
    # Load geometry externally
    geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

    # Run molecular dynamics
    md(geom, program='xtb', fmt='mol')


class TestXTBWrapper:

    def test_init(self, xtb):
        # Check class initialization
        assert isinstance(xtb, XTBWrapper)

        # Check temp directory created
        assert os.path.exists(xtb.temp_dir.name)

        # Clean up
        xtb.temp_dir.cleanup()

    @pytest.mark.parametrize('path,expected',
                             [('resources/geom_test.mol', 'geom_test')])
    def test_set_geometry(self, xtb, path, expected):
        # Load geometry externally
        geom = isicle.geometry.load(localfile(path))

        # Set geometry
        xtb.set_geometry(geom)

        # Check class instance
        assert isinstance(xtb.geom, Geometry)

        # Check basename attribute
        assert xtb.geom.basename == expected

        # Clean up
        xtb.temp_dir.cleanup()

    @pytest.mark.parametrize('fmt',
                             [('xyz'),
                              ('pdb'),
                              ('mol')])
    def test_save_geometry(self, xtb, fmt):
        # Load geometry externally
        geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

        # Set geometry
        xtb.set_geometry(geom)

        # Save geometry
        xtb.save_geometry(fmt=fmt)

        # Check if file exists
        assert os.path.exists(os.path.join(xtb.temp_dir.name,
                                           '{}.{}'.format(xtb.geom.basename,
                                                          fmt.lower())))

        # Clean up
        xtb.temp_dir.cleanup()

    # TODO: Get all combinations of tasks
    # TODO: meaningful assertions (currently just runs without error)
    @pytest.mark.parametrize('tasks,forcefield,optlevel',
                             [('optimize', 'gff', 'Normal'),
                              ('crest', 'gff', 'Normal'),
                              ('protonate', 'gff', 'Normal'),
                              (['optimize', 'crest'], ['gff', 'gfn2'], ['Normal', 'Tight']),
                              (['crest', 'protonate'], 'gfn2', 'Normal'),
                              (['optimize', 'protonate', 'deprotonate'], 'gff', 'Tight')])
    def test_job_type(self, xtb, tasks, forcefield, optlevel):
        # Load geometry externally
        geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

        # Set geometry
        xtb.set_geometry(geom)

        # Save geometry
        xtb.save_geometry(fmt='xyz')

        # Set up commandline
        xtb.job_type(tasks=tasks, forcefield=forcefield, optlevel=optlevel)

        # Clean up
        xtb.temp_dir.cleanup()

    def test_run(self, xtb):abns   
        # Load geometry externally
        geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

        # Set geometry
        xtb.set_geometry(geom)

        # Save geometry
        xtb.save_geometry(fmt='mol')

        # Set up commandline
        xtb.job_type()

        # Run
        xtb.run()

        # Clean up
        xtb.temp_dir.cleanup()

    def test_finish(self, xtb):
        # Load geometry externally
        geom = isicle.geometry.load(localfile('resources/geom_test.mol'))

        # Set geometry
        xtb.set_geometry(geom)

        # Copy example output to temp folder (i.e. assume nwchem was "run")
        shutil.copy2(localfile('resources/xtbopt.xyz'),
                     xtb.temp_dir.name)

        # Finish
        xtb.finish(keep_files=False)
