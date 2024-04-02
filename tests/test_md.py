import pytest
import isicle
from isicle.md import _program_selector, XTBWrapper, md, RDKitWrapper
from isicle.geometry import Geometry
import os
import shutil


def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)


@pytest.fixture()
def xtb():
    return XTBWrapper()


@pytest.fixture()
def rdw():
    return RDKitWrapper()


# test capitalization ignore
@pytest.mark.parametrize(
    "program,expected",
    [
        ("xtb", XTBWrapper),
        ("XTB", XTBWrapper),
        ("xTb", XTBWrapper),
        ("rdkit", RDKitWrapper),
        ("RDKit", RDKitWrapper),
        ("RDKIT", RDKitWrapper),
    ],
)
def test__program_selector(program, expected):
    # Check correct class is yielded
    assert isinstance(_program_selector(program), expected)


@pytest.mark.parametrize("program", [("xt"), ("amber"), ("rd")])
def test__program_selector_fail(program):
    # Check unsupported selections
    with pytest.raises(ValueError):
        _program_selector(program)


def test_md_xtb(xtb):
    # Load geometry externally
    geom = isicle.load(localfile("resources/geom_test_3D.mol"))

    # Run molecular dynamics
    md(
        geom,
        program="xtb",
        fmt="xyz",
        task="optimize",
        forcefield="gff",
        optlevel="Normal",
    )


def test_md_rdkit(rdw):
    # Load geometry externally
    geom = isicle.load(localfile("resources/geom_test_3D.mol"))

    # Run molecular dynamics
    md(
        geom,
        program="rdkit",
        method="etkdgv3",
    )


class TestXTBWrapper:
    def test_init(self, xtb):
        # Check class initialization
        assert isinstance(xtb, XTBWrapper)

        # Check temp directory created
        assert os.path.exists(xtb.temp_dir.name)

        # Clean up
        xtb.temp_dir.cleanup()

    @pytest.mark.parametrize(
        "path,expected", [("resources/geom_test_3D.mol", "geom_test_3D")]
    )
    def test_set_geometry(self, xtb, path, expected):
        # Load geometry externally
        geom = isicle.load(localfile(path))

        # Set geometry
        xtb.set_geometry(geom)

        # Check class instance
        assert isinstance(xtb.geom, Geometry)

        # Check basename attribute
        assert geom.basename == expected

        # Clean up
        xtb.temp_dir.cleanup()

    @pytest.mark.parametrize("fmt", [("xyz"), ("pdb"), ("mol")])
    def test_save_geometry(self, xtb, fmt):
        # Load geometry externally
        geom = isicle.load(localfile("resources/geom_test_3D.mol"))

        # Set geometry
        xtb.set_geometry(geom)

        # Save geometry
        xtb.save_geometry(fmt=fmt)

        # Check if file exists
        assert os.path.exists(
            os.path.join(
                xtb.temp_dir.name, "{}.{}".format(xtb.geom.basename, fmt.lower())
            )
        )

        # Clean up
        xtb.temp_dir.cleanup()

    @pytest.mark.parametrize(
        "tasks, forcefield, optlevel",
        [
            ("optimize", ["gff", "gfn2"], ["Normal", "Tight"]),
            (["crest", "protonate"], "gfn2", "Normal"),
        ],
    )
    def test_configure_failure(self, xtb, tasks, forcefield, optlevel):
        with pytest.raises(TypeError):
            xtb.configure(task=tasks, forcefield=forcefield, optlevel=optlevel)

    @pytest.mark.parametrize(
        "tasks,forcefield,optlevel",
        [
            ("optimize", "gff", "Normal"),
            ("crest", "gff", "Normal"),
            ("protonate", "gff", "Normal"),
        ],
    )
    def test_configure(self, xtb, tasks, forcefield, optlevel):
        # Load geometry externally
        geom = isicle.load(localfile("resources/geom_test_3D.mol"))

        # Set geometry
        xtb.set_geometry(geom)

        # Save geometry
        xtb.save_geometry(fmt="xyz")

        # Set up commandline
        xtb.configure(task=tasks, forcefield=forcefield, optlevel=optlevel)

        # Clean up
        xtb.temp_dir.cleanup()

    def test_run(self, xtb):
        # Load geometry externally
        geom = isicle.load(localfile("resources/geom_test_3D.mol"))

        # Set geometry
        xtb.set_geometry(geom)

        # Save geometry
        xtb.save_geometry(fmt="xyz")

        # Set up commandline
        xtb.configure()

        # Run
        xtb.run()

        # Clean up
        xtb.temp_dir.cleanup()

    def test_finish(self, xtb):
        # Load geometry externally
        geom = isicle.load(localfile("resources/geom_test_3D.mol"))

        # Set geometry
        xtb.set_geometry(geom)

        # Save geometry
        xtb.save_geometry(fmt="xyz")

        # Set up commandline
        xtb.configure()

        # Run
        xtb.run()

        # Finish
        xtb.finish()

        # Clean up
        xtb.temp_dir.cleanup()


class TestRDKitWrapper:
    def test_init(self, rdw):
        # Check class initialization
        assert isinstance(rdw, RDKitWrapper)

    @pytest.mark.parametrize(
        "path,expected", [("resources/geom_test_3D.mol", "geom_test_3D")]
    )
    def test_set_geometry(self, rdw, path, expected):
        # Load geometry externally
        geom = isicle.load(localfile(path))

        # Set geometry
        rdw.set_geometry(geom)

        # Check class instance
        assert isinstance(rdw.geom, Geometry)

        # Check basename attribute
        assert geom.basename == expected

    @pytest.mark.parametrize(
        "method, numConfs",
        [("distance", "test"), ("etkdg", "test")],
    )
    def test_configure_failure(self, rdw, method, numConfs):
        with pytest.raises(ValueError):
            rdw.configure(method=method, numConfs=numConfs)

    @pytest.mark.parametrize(
        "method, numConfs, pruneRmsThresh, forceTol, randomSeed",
        [
            ("distance", 10, -1.0, 0.001, -1),
        ],
    )
    def test_configure_distance(
        self, rdw, method, numConfs, pruneRmsThresh, forceTol, randomSeed
    ):
        # Load geometry externally
        geom = isicle.load(localfile("resources/geom_test_3D.mol"))

        # Set geometry
        rdw.set_geometry(geom)

        # Set up commandline
        rdw.configure(
            method=method,
            numConfs=numConfs,
            pruneRmsThresh=pruneRmsThresh,
            forceTol=forceTol,
            randomSeed=randomSeed,
        )

    @pytest.mark.parametrize(
        "method, numConfs",
        [("etkdg", 10), ("etdg", 10)],
    )
    def test_configure_et(self, rdw, method, numConfs):
        # Load geometry externally
        geom = isicle.load(localfile("resources/geom_test_3D.mol"))

        # Set geometry
        rdw.set_geometry(geom)

        # Set up commandline
        rdw.configure(method=method, numConfs=numConfs)

    def test_submit(self, rdw):
        # Load geometry externally
        geom = isicle.load(localfile("resources/geom_test_3D.mol"))

        # Set geometry
        rdw.set_geometry(geom)

        # Set up commandline
        rdw.configure()

        # Run
        rdw.submit()

    def test_finish(self, rdw):
        # Load geometry externally
        geom = isicle.load(localfile("resources/geom_test_3D.mol"))

        # Set geometry
        rdw.set_geometry(geom)

        # Set up commandline
        rdw.configure()

        # Run
        rdw.submit()

        # Finish
        rdw.finish()
