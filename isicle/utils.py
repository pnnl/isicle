import numpy as np
import pandas as pd


def safelist(x):
    """
    Ensures passed object is of correct list-like format.

    Parameters
    ----------
    x : any
        Object to be cast as list.
    Returns
    -------
    out : list_like
        Input safely cast to list-like.

    """

    if not isinstance(x, (list, pd.core.series.Series, np.ndarray)):
        return [x].copy()
    return x.copy()
