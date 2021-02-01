import pytest
import isicle


@pytest.fixture()
def nwc():
    return isicle.qm.NWChemWrapper()


class TestNWChemWrapper:

    def test_init(self, nwc):
        raise NotImplementedError

    # Note: mostly relies on `geometry.load`, so
    # no need to test fully. Just make sure attributes
    # are assigned and outputs generated correctly.
    def test_load_geometry(self, nwc):
        raise NotImplementedError

    def test__set_template(self, nwc):
        raise NotImplementedError

    def test_configure(self, nwc):
        raise NotImplementedError

    # Not particularly testable
    def test_run(self, nwc):
        raise NotImplementedError

    def test_finish(self, nwc):
        raise NotImplementedError