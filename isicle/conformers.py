from statsmodels.stats.weightstats import DescrStatsW
import pandas as pd
from isicle.geometry import Geometry, MDOptimizedGeometry, DFTOptimizedGeometry


def _method_selector(method):
    method_map = {'boltzmann': boltzmann,
                  'simple': simple_average,
                  'lowest': lowest_energy,
                  'threshold': threshold}

    if method.lower() in method_map.keys():
        return method_map[method.lower()]
    else:
        raise ValueError('{} not a supported conformer reduction method.'.format(method))


def reduce(value, method='boltzmann', **kwargs):
    f = _method_selctor(method)
    return f(value, energy, **kwargs)


def boltzmann(value, energy=None, index=None):
    if index is None:
        index = -1

    df = pd.DataFrame({'value': value, 'energy': energy, 'index': index})

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
        return res.drop(columns='index')

    return res


def simple_average(value, index=None):
    if index is None:
        index = -1

    df = pd.DataFrame({'value': value, 'index': index})

    res = df.groupby(['index'], as_index=False).mean()

    if index is None:
        return res.drop(columns='index')

    return res


def lowest_energy(value, energy=None, index=None):
    if index is None:
        index = -1

    df = pd.DataFrame({'value': value, 'energy': energy, 'index': index})

    res = df.loc[df.groupby('index')['energy'].idxmin()]

    if index is None:
        return res.drop(columns='index')

    return res


def threshold(value, energy=None, threshold=5, index=None):
    if index is None:
        index = -1

    df = pd.DataFrame({'value': value, 'energy': energy, 'index': index})

    df = df.loc[df['energy'] <= threshold, :]

    res = df.groupby(['index'], as_index=False).mean()

    if index is None:
        return res.drop(columns='index')

    return res


def _are_Geometry_instances(objects):
    return all(isinstance(x, (Geometry, MDOptimizedGeometry, DFTOptimizedGeometry)) for x in objects)


def build_conformational_ensemble(geometries):
    if _are_Geometry_instances(geometries) is True:
        return ConformationalEnsemble(geometries)
    else:
        raise TypeError('Conformers must be of type Geometry or related subclass.')


class ConformationalEnsemble(list):
    def __init__(self):
        pass

    def reduce(self):
        raise NotImplementedError()

    def _apply_method(self, method, **kwargs):
        # Check for attribute
        if not all(hasattr(x, method) for x in self):
            raise AttributeError('{} not found for all conformational sample members.')

        # Apply method to collection
        result = [getattr(x, method)(**kwargs) for x in self]

        # Return ConformationalEnsemble if correct result type
        if _are_Geometry_instances(result) is True:
            return ConformationalEnsemble(result)

        # Return result as-is
        else:
            return result

    def _apply_function(self, func, **kwargs):
        # Apply method to collection
        result = [func(x, **kwargs) for x in self]

        # Return ConformationalEnsemble if correct result type
        if _are_Geometry_instances(result) is True:
            return ConformationalEnsemble(result)

        # Return result as-is
        else:
            return result

    def apply(self, func=None, method=None, **kwargs):
        if func is not None:
            return self._apply_function(func, **kwargs)

        if method is not None:
            return self._apply_method(method, **kwargs)

        raise ValueError('Must supply `func` or `method`.')
