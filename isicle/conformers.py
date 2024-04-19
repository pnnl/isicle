import numpy as np
import pandas as pd
from statsmodels.stats.weightstats import DescrStatsW
from os.path import join

from isicle import io
from isicle.geometry import Geometry, XYZGeometry
from isicle.utils import TypedList, safelist, mkdtemp
from rdkit import Chem
from rdkit.Chem import PropertyMol
from confpass import confpass


def _function_selector(func):
    """
    Selects a supported reduction function for reducing a set of conformers.

    Parameters
    ----------
    func : str
        Alias for function selection (one of "boltzmann", "simple", "lowest",
        or "threshold").

    Returns
    -------
    func
        Conformer reduction function.
    """
    # Mapping between names and functions
    func_map = {
        "boltzmann": boltzmann,
        "simple": simple_average,
        "lowest": lowest_energy,
        "threshold": energy_threshold,
    }

    # Check for function by name
    if func.lower() in func_map:
        return func_map[func.lower()]
    else:
        # Not a named/implemented function
        raise ValueError("{} not a supported reduction function.".format(func))


def _energy_based(func):
    """
    Checks whether function employs an energy-based reduction operation.

    Parameters
    ----------
    func : function
        Conformer reduction function.

    Returns
    -------
    bool
        True if energy based, otherwise False.
    """
    # Check if among energy-based functions
    if func in [boltzmann, lowest_energy, energy_threshold]:
        return True

    # Not energy based
    return False


def reduce(value, func="boltzmann", **kwargs):
    """
    Combine values according to indicated function.

    Parameters
    ----------
    value : :obj:`~numpy.array`
        Array containing values that will be combined.
    func : str
        Alias for function selection (one of "boltzmann", "simple", "lowest",
        or "threshold").
    kwargs
        Additional keyword arguments passed to `func`.

    Returns
    -------
    :obj:`~pandas.DataFrame`
        Result of reduction operation.

    """
    # Select function
    f = _function_selector(func)

    # Energy-based method
    if _energy_based(f):
        energy = kwargs.pop("energy")
        return f(value, energy, **kwargs)

    # Other method
    return f(value, **kwargs)


def boltzmann(value, energy, index=None, atom=None):
    """
    Combine values according to a Boltzmann-weighted average.

    Parameters
    ----------
    value : :obj:`~numpy.array`
        Array containing values that will be combined.
    energy : :obj:`~numpy.array`
        Array containing energy values that correspond to entries in `value`.
    index : None or :obj:`~numpy.array`
        Index by which to group values for averaging.
    atom : None or :obj:`~numpy.array`
        Atom by which to group values for averaging.

    Returns
    -------
    :obj:`~pandas.DataFrame`
        Result of reduction operation.
    """

    # Placeholder for index
    if index is None:
        index = np.full_like(value, -1)

    # Placeholder for atom
    if atom is None:
        atom = np.full_like(value, -1)

    # Initialize data frame
    df = pd.DataFrame.from_dict(
        {"value": value, "energy": energy, "index": index, "atom": atom}
    )

    # Result container
    res = []

    # Iterate over unique indices
    for name, group in df.groupby(["index", "atom"]):
        # Compute relative delta G
        g = group["energy"] * 627.503
        mn = g.min()
        relG = g - mn

        # Compute Boltzmann weighting factors
        b = np.exp(-relG / 0.5924847535)
        w = (b / b.sum()) * len(b)

        # Compute weighted statistics
        ws = DescrStatsW(group["value"], weights=w, ddof=0)

        # Append to container
        res.append([name[0], name[1], ws.mean, ws.std, len(group.index)])

    # Initialize data frame
    res = pd.DataFrame(res, columns=["index", "atom", "mean", "std", "n"])

    # Drop index if not supplied
    if np.all(index == -1):
        return res.drop(columns=["index", "atom"]).iloc[0]

    return res


def simple_average(value, index=None, atom=None):
    """
    Combine values according to a simple average.

    Parameters
    ----------
    value : :obj:`~numpy.array`
        Array containing values that will be combined.
    index : None or :obj:`~numpy.array`
        Index by which to group values for averaging.
    atom : None or :obj:`~numpy.array`
        Atom by which to group values for averaging.

    Returns
    -------
    :obj:`~pandas.DataFrame`
        Result of reduction operation.
    """

    # Placeholder for index
    if index is None:
        index = np.full_like(value, -1)

    # Placeholder for atom
    if atom is None:
        atom = np.full_like(value, -1)

    # Initialize data frame
    df = pd.DataFrame.from_dict({"value": value, "index": index, "atom": atom})

    # Average per unique index
    res = df.groupby(["index", "atom"], as_index=False).agg(
        {"value": ["mean", "std", "count"]}
    )
    # Rename columns
    res.columns = ["index", "atom", "mean", "std", "n"]

    # Drop indices if not supplied
    if np.all(index == -1):
        return res.drop(columns=["index", "atom"]).iloc[0]

    return res


def lowest_energy(value, energy, index=None, atom=None):
    """
    Combine values according to lowest energy.

    Parameters
    ----------
    value : :obj:`~numpy.array`
        Array containing values that will be combined.
    energy : :obj:`~numpy.array`
        Array containing energy values that correspond to entries in `value`.
    index : None or :obj:`~numpy.array`
        Index by which to group values for averaging.
    atom : None or :obj:`~numpy.array`
        Atom by which to group values for averaging.

    Returns
    -------
    :obj:`~pandas.DataFrame`
        Result of reduction operation.
    """

    # Placeholder for index
    if index is None:
        index = np.full_like(value, -1)

    # Placeholder for atom
    if atom is None:
        atom = np.full_like(value, -1)

    # Initialize data frame
    df = pd.DataFrame.from_dict(
        {"value": value, "energy": energy, "index": index, "atom": atom}
    )
    # Take minimum energy per unique index
    res = df.loc[df.groupby(["index", "atom"])["energy"].idxmin()]

    # Drop indices if not supplied
    if np.all(index == -1):
        return res.drop(columns=["index", "atom"]).iloc[0]

    return res


def energy_threshold(value, energy, threshold=5, index=None, atom=None):
    """
    Combine values with energy below a given threshold according to a simple
    average.

    Parameters
    ----------
    value : :obj:`~numpy.array`
        Array containing values that will be combined.
    energy : :obj:`~numpy.array`
        Array containing energy values that correspond to entries in `value`.
    index : None or :obj:`~numpy.array`
        Index by which to group values for averaging.
    atom : None or :obj:`~numpy.array`
        Atom by which to group values for averaging.

    Returns
    -------
    :obj:`~pandas.DataFrame`
        Result of reduction operation.
    """

    # Placeholder for index
    if index is None:
        index = np.full_like(value, -1)

    # Placeholder for atom
    if atom is None:
        atom = np.full_like(value, -1)

    # Initialize data frame
    df = pd.DataFrame.from_dict(
        {"value": value, "energy": energy, "index": index, "atom": atom}
    )

    # Filter by energy
    df = df.loc[df["energy"] <= threshold, :]

    # Aggregate
    res = df.groupby(["index", "atom"], as_index=False).agg(
        {"value": ["mean", "std", "count"]}
    )

    # Rename columns
    res.columns = ["index", "atom", "mean", "std", "n"]

    # Drop indices if not supplied
    if index is None:
        return res.drop(columns=["index", "atom"]).iloc[0]

    return res


def transform(
    value, m={"H": 1.0, "C": 1.0}, b={"H": 0.0, "C": 0.0}, index=None, atom=None
):
    """
    Perform linear transformation with values using provided parameters.

    Parameters
    ----------
    value : :obj: `~numpy.array`
        Array containing vales that will be transformed.
    m : float or dict
        Slope value
    b : float or dict
        Y-intercept value
    index : None or :obj: `~numpy.array`
        Index by which to group values for transforming.
    atom : None or :obj: `~numpy.array`
        Atom by which to group values for transforming.

    Returns
    -------
    :obj: `~pandas.DataFrame`
        Result of transformation operation.
    """

    # Placeholder for index
    if index is None:
        index = np.full_like(value, -1)

    # Placeholder for atom
    if atom is None:
        atom = np.full_like(value, -1)

    # Initialize data frame
    df = pd.DataFrame.from_dict({"value": value, "index": index, "atom": atom})

    # Process with per-atom values
    if isinstance(m, dict):
        res = pd.DataFrame()
        for idx in m:
            part = df.loc[df["atom"] == idx].copy()
            part["new_value"] = part["value"].apply(lambda x: m[idx] * x + b[idx])

            res = pd.concat([res, part])

    # Process with global values
    else:
        res = df.copy()
        res["new_value"] = res["value"].apply(lambda x: m * x + b)

    return res


def build_conformational_ensemble(geometries):
    """
    Create a conformational ensemble from a collection of geometries.

    Parameters
    ----------
    geometries : list of :obj:`~isicle.geometry.Geometry` or related subclass
        Collection of geometry instances.

    Returns
    -------
    :obj:`~isicle.conformers.ConformationalEnsemble`
        Conformational ensemble.
    """

    return ConformationalEnsemble(geometries)


class ConformationalEnsemble(TypedList):
    """
    Collection of :obj:`~isicle.geometry.Geometry`, or related subclass,
    instances.
    """

    def __init__(self, *args):
        """
        Initialize :obj:`~isicle.conformers.ConformationalEnsemble` instance.

        Parameters
        ----------
        *args
            Objects to comprise the conformational ensemble.

        """

        super().__init__((Geometry, XYZGeometry), *args)

    def _check_attributes(self, attr):
        """
        Check if all ensemble members have the supplied attribute.

        Parameters
        ----------
        attr : str
            Attribute to check.

        Raises
        ------
        AttributeError
            If all members do not have `attr`.

        """

        value = [x.get___dict__() for x in self]
        for key in safelist(attr):
            if not all(key in x for x in value):
                raise AttributeError(
                    '"{}" not found for all conformational '
                    "ensemble members.".format(attr)
                )
            value = [x.get(key) for x in value]

        # if type(value[0]) is dict:
        #     raise AttributeError('"{}" has additional keys: {}'.format(attr, value[0].keys()))

    def reduce(self, attr, func="boltzmann", **kwargs):
        """
        Combine attribute values according to indicated function.

        Parameters
        ----------
        attr : str
            Attribute that will be combined.
        func : str
            Alias for function selection (one of "boltzmann", "simple",
            "lowest", or "threshold").
        kwargs
            Additional keyword arguments passed to `func`.

        Returns
        -------
        :obj:`~pandas.DataFrame`
            Result of reduction operation.

        """

        # Select reduction function
        f = _function_selector(func)

        # Check for primary attribute
        self._check_attributes(attr)

        # Check for energy attribute
        if _energy_based(f):
            self._check_attributes("energy")

        # Extract (possibly nested) value attribute
        value = [x.get___dict__() for x in self]
        for key in safelist(attr):
            value = [x.get(key) for x in value]

        # Check nested values
        if isinstance(value[0], dict):
            # Check index
            if "index" in value[0]:
                index = np.array([x["index"] for x in value]).flatten()
                pad = int(len(index) / len(self))
            else:
                index = None
                pad = 1

            # Check atom
            if "atom" in value[0]:
                atom = np.array([x["atom"] for x in value]).flatten()
            else:
                atom = None

            # Special case for CCS
            if "mean" in value[0] and "std" in value[0]:
                value = np.array([x["mean"] for x in value]).flatten()

            else:
                value = np.array([x[attr] for x in value]).flatten()

        # Not nested
        else:
            index = None
            atom = None
            pad = 1

        # Extract energy attribute
        if _energy_based(f):
            energy = np.array(
                [np.repeat(x.get___dict__()["energy"], pad) for x in self]
            )
            energy = energy.flatten()

            # Exectue energy-based method
            return f(value, energy, index=index, atom=atom, **kwargs)

        # Execute other method
        return f(value, index=index, atom=atom, **kwargs)

    def _apply_method(self, method, **kwargs):
        """
        Process conformational ensemble members according to supplied method.

        Parameters
        ----------
        method : str
            Method by which ensemble members will be processed.
        kwargs
            Keyword arguments passed to `method`.

        Returns
        -------
        :obj:`~isicle.conformers.ConformationalEnsemble` or list
            Result of operation, type depends on `method` return type.

        """

        # Check for attribute
        if not all(hasattr(x, method) for x in self):
            raise AttributeError(
                '"{}" not found for all conformational '
                "ensemble members.".format(method)
            )

        # Apply method to collection
        result = [getattr(x, method)(**kwargs) for x in self]

        # Return ConformationalEnsemble if correct result type
        try:
            return ConformationalEnsemble(result)

        # Return result as-is
        except:
            return result

    def _apply_function(self, func, **kwargs):
        """
        Process conformational ensemble members according to supplied function.

        Parameters
        ----------
        func : function
            Function by which ensemble members will be processed.
        kwargs
            Keyword arguments passed to `func`.

        Returns
        -------
        :obj:`~isicle.conformers.ConformationalEnsemble` or list
            Result of operation, type depends on `func` return type.

        """

        # Apply method to collection
        result = [func(x, **kwargs) for x in self]

        # Return ConformationalEnsemble if correct result type
        try:
            return ConformationalEnsemble(result)

        # Return result as-is
        except:
            return result

    def apply(self, func=None, method=None, **kwargs):
        """
        Process conformational ensemble members according to supplied function
        or method.

        Parameters
        ----------
        func : function
            Function by which ensemble members will be processed.
        method : str
            Method by which ensemble members will be processed.
        kwargs
            Keyword arguments passed to `method`.

        Returns
        -------
        :obj:`~isicle.conformers.ConformationalEnsemble` or list
            Result of operation, type depends on `method` return type.

        Raises
        ------
        ValueError
            If neither `func` nor `method` is supplied.

        """

        # Apply function
        if func is not None:
            return self._apply_function(func, **kwargs)

        # Apply method
        if method is not None:
            return self._apply_method(method, **kwargs)

        raise ValueError("Must supply `func` or `method`.")

    def get_structures(self):
        """
        Extract all structures from containing object as a conformational ensemble.

        Returns
        -------
        :obj:`~isicle.conformers.ConformationalEnsemble`
            Conformational ensemble.

        """

        # Check for geom attribute
        self._check_attributes("geom")

        # Build and return
        return build_conformational_ensemble([x.geom for x in self])


class ConfpassEnsemble(ConformationalEnsemble):
    """
    Collection of :obj:`~isicle.geometry.Geometry`, or related subclass,
    instances. This class REQUIRES the software CONFPASS to be installed
    in Python path.
    """

    _defaults = ["temp_dir", "basename", "sdf_file", "priority", "cp"]
    _default_value = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args)
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))
        self.__dict__.update(**kwargs)

    def _create_PropertyMolID(self, geom):
        """
        Create molecule identifiers from attributes.
        """
        PropertyMolID = []
        try:
            # geom._check_attributes('basename')
            PropertyMolID.append(geom.basename)
        except AttributeError:
            pass
        try:
            # geom._check_attributes('adduct')
            PropertyMolID.append(geom.adduct)
        except AttributeError:
            pass
        try:
            # geom._check_attributes('adductID')
            PropertyMolID.append(geom.adductID)
        except AttributeError:
            pass
        try:
            # geom._check_attributes('conformerID')
            PropertyMolID.append(geom.conformerID)
        except AttributeError:
            pass
        PropertyMolID = [str(i) for i in PropertyMolID]
        return "_".join(PropertyMolID)

    def _make_geometry(self):
        """
        Create RDKit PropertyMol representations of molecule.
        """
        pmObjs = [PropertyMol.PropertyMol(geom.mol) for geom in self]
        pmIDs = [self._create_PropertyMolID(geom) for geom in self]
        for pm, pmID in zip(pmObjs, pmIDs):
            pm.SetProp("_Name", pmID)
        return pmObjs

    def _set_basename(self):
        """
        Ensure and set the basename is same for all members of conformation ensemble.
        """
        self._check_attributes("basename")
        value = list(set([x.basename for x in self]))
        assert len(value) == 1, "all basenames are not equal"
        self.basename = value[0]

    def save_geometry(self):
        """
        Save internal conformational ensemble geometry representation to SDF file.

        """
        # get basename for conformational ensemble
        self._set_basename()

        # create PropertyMol RDKit versions
        pms = self._make_geometry()
        # create, save temporary multi-molecule sdf file
        self.temp_dir = mkdtemp()
        self.sdf_file = join(self.temp_dir, "{}.{}".format(self.basename, "sdf"))
        io.save(self.sdf_file, pms)

    def _apply_preopt(self, method, x=0.8, x_as=0.2, n=3):
        """
        Run the pre-optimization CONFPASS module.
        Parameters
        ----------
        method : str
            Priority list assembling method
        x : float
            Hyperparameter for the pipe_x method, default=0.8.
            Select clustering result when n_clusters=total_conformer_no*x
        x_as : float
            Hyperparameter for the pipe_x_as method, default=0.2.
            Ex: 20% (pipe x priority list) + 80% (pipe d priority list)
        n : int
            Hyperparamter for the nth method, default = 3.
            Prioritise every nth conformer.
        """
        cp = confpass.conp([self.sdf_file])
        cp.get_priority(method=method, x=x, x_as=x_as, n=n)
        self.cp = cp
        self.priority = cp.priority_df.loc[0, "priority_ls"]

    def _apply_reopt(self, method, x=0.8, x_as=0.2, n=3):
        """
        Run the post-optimization CONFPASS module.
        """
        cp = confpass.pas([self.sdf_file])
        cp.get_priority(method=method, x=x, x_as=x_as, n=n)
        self.cp = cp

    def _check_module(self, module):
        """
        Check specified CONFPASS module is supported
        """
        module = module.lower()
        assert module in ["conp", "pas"], "Specified CONFPASS module not supported"
        return module

    def _check_method(self, method):
        """
        Check specified CONFPASS method is supported
        """
        method = method.lower()
        assert method in [
            "pipe_x_as",
            "pipe_x",
            "pipe_as",
            "nth",
            "random",
            "ascend",
        ], "Specified CONFPASS method not supported"
        return method

    def apply(self, module, method, **kwargs):
        """
        Parameters
        ----------
        module : str
            Specify module as
            `conp` pre-reoptimisations: clustering and priority list generation
            `pas` post-reoptimisations: predict the completeness of the reoptimisation process
        method : str
            `pipe_x_as`, `pipe_x`,`pipe_as`,`nth`,`random`,`ascend`

        """
        module = self._check_module(module)
        method = self._check_method(method)

        if module == "conp":
            self._apply_preopt(method, **kwargs)
        elif module == "pas":
            self._apply_reopt(method, **kwargs)
