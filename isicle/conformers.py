from statsmodels.stats.weightstats import DescrStatsW
import pandas as pd
from isicle.geometry import Geometry, MDOptimizedGeometry, DFTOptimizedGeometry
import numpy as np
from isicle.utils import TypedList


def _function_selector(func):
    '''
    Selects a supported reduction function for reducing a set of conformers.

    Parameters
    ----------
    function : str
        Alias for function selection (one of "boltzmann", "simple", "lowest",
        or "threshold").

    Returns
    -------
    func
        Conformer reduction function.

    '''

    func_map = {'boltzmann': boltzmann,
                'simple': simple_average,
                'lowest': lowest_energy,
                'threshold': threshold}

    if func.lower() in func_map.keys():
        return func_map[func.lower()]
    else:
        raise ValueError('{} not a supported reduction function.'.format(func))


def _energy_based(f):
    if f in [boltzmann, lowest_energy, threshold]:
        return True

    return False


def reduce(value, func='boltzmann', **kwargs):
    f = _function_selector(func)

    # Energy-based method
    if _energy_based(f):
        energy = kwargs.pop('energy')
        return f(value, energy, **kwargs)

    # Other method
    return f(value, **kwargs)


def boltzmann(value, energy, index=None):
    df = pd.DataFrame({'value': value, 'energy': energy, 'index': -1})

    if index is not None:
        df['index'] = index

    res = []
    for name, group in df.groupby(['index']):
        g = group['energy'] * 627.503
        mn = g.min()
        relG = g - mn
        b = np.exp(-relG / 0.5924847535)
        w = (b / b.sum()) * len(b)

        ws = DescrStatsW(group['value'], weights=w, ddof=0)

        res.append([name, ws.mean, ws.std, len(group.index)])

    res = pd.DataFrame(res, columns=['index', 'mean', 'std', 'n'])

    if index is None:
        return res.drop(columns='index').iloc[0]

    return res


def simple_average(value, index=None):
    df = pd.DataFrame({'value': value, 'index': -1})

    if index is not None:
        df['index'] = index

    res = df.groupby(['index'], as_index=False).agg({'value':
                                                     ['mean', 'std', 'count']})
    res.columns = ['index', 'mean', 'std', 'n']

    if index is None:
        return res.drop(columns='index').iloc[0]

    return res


def lowest_energy(value, energy, index=None):
    df = pd.DataFrame({'value': value, 'energy': energy, 'index': -1})

    if index is not None:
        df['index'] = index

    res = df.loc[df.groupby('index')['energy'].idxmin()]

    if index is None:
        return res.drop(columns='index').iloc[0]

    return res


def threshold(value, energy, threshold=5, index=None):
    df = pd.DataFrame({'value': value, 'energy': energy, 'index': -1})

    if index is not None:
        df['index'] = index

    df = df.loc[df['energy'] <= threshold, :]

    res = df.groupby(['index'], as_index=False).agg({'value':
                                                     ['mean', 'std', 'count']})
    res.columns = ['index', 'mean', 'std', 'n']

    if index is None:
        return res.drop(columns='index').iloc[0]

    return res


def build_conformational_ensemble(geometries):
    return ConformationalEnsemble(geometries)


class ConformationalEnsemble(TypedList):
    def __init__(self, *args):
        super().__init__((Geometry, MDOptimizedGeometry, DFTOptimizedGeometry),
                         *args)

    def _check_attributes(self, attr):
        if not all(hasattr(x, attr) for x in self):
            raise AttributeError('"{}" not found for entire conformational'
                                 'sample members.'.format(attr))

    def reduce(self, attr, func='boltzmann', index=False, **kwargs):
        f = _function_selector(func)

        # Check for primary attribute
        self._check_attributes(attr)

        # Check for energy attribute
        if _energy_based(f):
            self._check_attributes('energy')

        # Check for index attribute
        if index is True:
            self._check_attributes('index')

            # Extract attribute
            index = np.array([getattr(x, 'index') for x in self]).flatten()
            pad = int(len(index) / len(self))

        # No index
        else:
            index = None
            pad = 1

        # Extract value attribute
        value = np.array([getattr(x, attr) for x in self]).flatten()

        # Extract energy attribute
        if _energy_based(f):
            energy = np.array([[getattr(x, 'energy')] * pad for x in self])
            energy = energy.flatten()

            # Exectue energy-based method
            return f(value, energy, index=index, **kwargs)

        # Execute other method
        return f(value, index=index, **kwargs)

    def _apply_method(self, method, **kwargs):
        # Check for attribute
        self._check_attributes(method)

        # Apply method to collection
        result = [getattr(x, method)(**kwargs) for x in self]

        # Return ConformationalEnsemble if correct result type
        try:
            return ConformationalEnsemble(result)

        # Return result as-is
        except:
            return result

    def _apply_function(self, func, **kwargs):
        # Apply method to collection
        result = [func(x, **kwargs) for x in self]

        # Return ConformationalEnsemble if correct result type
        try:
            return ConformationalEnsemble(result)

        # Return result as-is
        except:
            return result

    def apply(self, func=None, method=None, **kwargs):
        if func is not None:
            return self._apply_function(func, **kwargs)

        if method is not None:
            return self._apply_method(method, **kwargs)

        raise ValueError('Must supply `func` or `method`.')
