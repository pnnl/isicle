import numpy as np
import pandas as pd
import collections
from pkg_resources import resource_filename


def safelist(x):
    '''
    Ensures passed object is of correct format.

    Parameters
    ----------
    x : any
        Object to be cast as list.
    Returns
    -------
    list, :obj:`~pd.core.series.Series`, or :obj:`~np.ndarray`
        Input safely cast to list-like.

    '''

    if not isinstance(x, (list, pd.core.series.Series, np.ndarray)):
        return [x].copy()
    return x.copy()


class TypedList(collections.abc.MutableSequence):
    '''
    Mutable sequence that requires all members be of select type(s).

    Attributes
    ----------
    oktypes : type or list of types
        Object types allowed list membership.
    list : list
        Internal list representation.

    '''

    def __init__(self, oktypes, *args):
        '''
        Initialize :obj:`~isicle.utils.TypedList` instance.

        Parameters
        ----------
        oktypes : type or list of types
            Object types allowed list membership.
        *args
            Objects to comprise the type-restricted list.

        '''

        self.oktypes = oktypes
        self.list = list()
        try:
            self.extend(*args)
        except:
            self.extend(list(args))

    def check(self, v):
        '''
        Check if supplied value is of allowed type(s).

        Raises
        ------
        TypeError
            If value is not of allowed type(s).

        '''

        if not isinstance(v, self.oktypes):
            raise TypeError(v)

    def __len__(self):
        return len(self.list)

    def __getitem__(self, i):
        return self.list[i]

    def __delitem__(self, i):
        del self.list[i]

    def __setitem__(self, i, v):
        self.check(v)
        self.list[i] = v

    def insert(self, i, v):
        self.check(v)
        self.list.insert(i, v)

    def __str__(self):
        return str(self.list)

    def __repr__(self):
        return self.__str__()


def atomic_masses():
    path = resource_filename('isicle', 'resources/atomic_masses.tsv')
    return pd.read_csv(path, sep='\t')
