import os


def localfile(path):
    '''
    Returns path relative to this file.
    Parameters
    ----------
    path : str
        Relative path to file.
    Returns
    -------
    path : str
        Absolute path.
    '''

    return os.path.join(os.path.dirname(__file__), path)
