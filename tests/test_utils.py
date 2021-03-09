import pytest
import isicle
import numpy as np
from pandas.core.series import Series


@pytest.fixture()
def tlint():
    return isicle.utils.TypedList(int, 1, 2, 3, 4)


@pytest.mark.parametrize('x,expected',
                         [('a', ['a']),
                          (['a', 'b', 'c'], ['a', 'b', 'c']),
                          (1, [1]),
                          ([1, 2, 3], [1, 2, 3])])
def test_safelist(x, expected):
    # list
    assert isicle.utils.safelist(x) == expected

    # array
    assert np.all(isicle.utils.safelist(np.array(x)) == np.array(expected))

    # series
    assert (Series(x) == Series(expected)).all()


class TestTypedList:
    @pytest.mark.parametrize('allowed,values',
                             [(int, [1, 2, 3, 4]),
                              (str, ['1', '2', '3', '4']),
                              ((int, str), [1, '2', 3, '4'])])
    def test_init(self, allowed, values):
        isicle.utils.TypedList(allowed, *values)
        isicle.utils.TypedList(allowed, values)

    @pytest.mark.parametrize('allowed,values',
                             [(int, ['1', 2, 3, 4]),
                              (str, [1, '2', '3', '4']),
                              ((int, str), [1, '2', 3, 1.1])])
    def test_init_fail(self, allowed, values):
        with pytest.raises(TypeError):
            isicle.utils.TypedList(allowed, *values)
        with pytest.raises(TypeError):
            isicle.utils.TypedList(allowed, values)

    @pytest.mark.parametrize('value',
                             [(2)])
    def test_check(self, tlint, value):
        tlint.check(value)

    @pytest.mark.parametrize('value',
                             [('2')])
    def test_check_fail(self, tlint, value):
        with pytest.raises(TypeError):
            tlint.check(value)

    @pytest.mark.parametrize('allowed,values,expected',
                             [(int, [1, 2, 3, 4], 4),
                              (str, ['1', '2', '3', '4'], 4),
                              ((int, str), [1, '2', 3, '4'], 4)])
    def test_len(self, allowed, values, expected):
        assert len(isicle.utils.TypedList(allowed, values)) == expected

    @pytest.mark.parametrize('index,expected',
                             [(2, 3),
                              (3, 4),
                              (0, 1)])
    def test_getitem(self, tlint, index, expected):
        assert tlint[index] == expected

    @pytest.mark.parametrize('index,expected',
                             [(2, [1, 2, 4]),
                              (3, [1, 2, 3]),
                              (0, [2, 3, 4])])
    def test_delitem(self, tlint, index, expected):
        del tlint[index]
        assert tlint.list == expected

    @pytest.mark.parametrize('index,value,expected',
                             [(0, 9, [9, 2, 3, 4]),
                              (1, 9, [1, 9, 3, 4]),
                              (-1, 9, [1, 2, 3, 9])])
    def test_setitem(self, tlint, index, value, expected):
        tlint[index] = value
        assert tlint.list == expected

    @pytest.mark.parametrize('index,value,expected',
                             [(0, 9, [9, 1, 2, 3, 4]),
                              (4, 9, [1, 2, 3, 4, 9]),
                              (-1, 9, [1, 2, 3, 9, 4])])
    def test_insert(self, tlint, index, value, expected):
        tlint.insert(index, value)
        assert tlint.list == expected

    @pytest.mark.parametrize('expected',
                             [('[1, 2, 3, 4]')])
    def test_str(self, tlint, expected):
        assert str(tlint) == expected

    @pytest.mark.parametrize('expected',
                             [('[1, 2, 3, 4]')])
    def test_repr(self, tlint, expected):
        assert repr(tlint) == expected
