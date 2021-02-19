import numpy as np
import pandas as pd


def safelist(x):
    """
    Ensures passed object is of correct format.

    Parameters
    ----------
    x : any
        Object to be cast as list.
    Returns
    -------
    list, :obj:`~pd.core.series.Series`, or :obj:`~np.ndarray`
        Input safely cast to list-like.

    """

    if not isinstance(x, (list, pd.core.series.Series, np.ndarray)):
        return [x].copy()
    return x.copy()
