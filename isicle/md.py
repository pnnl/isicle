import glob
import os
import subprocess
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import rdDistGeom

import isicle
from isicle.geometry import Geometry
from isicle.interfaces import WrapperInterface

"""
Files resulting from an xtb job always run in the same directory that the command is
issued in, no matter where the input is. Can direct the .log file, but no other files.
"""


def _backend_selector(backend):
    """
    Selects a supported molecular dynamics backend for associated simulation.
    Currently only NWChem has been implemented.

    Parameters
    ----------
    backend : str
        Alias for backend selection (xtb).

    Returns
    -------
    backend
        Wrapped functionality of the selected backend. Must implement
        :class:`~isicle.interfaces.WrapperInterface`.

    """

    backend_map = {"xtb": XTBWrapper, "rdkit": RDKitWrapper}

    if backend.lower() in backend_map.keys():
        return backend_map[backend.lower()]()
    else:
        raise ValueError(
            "{} not a supported molecular dynamics backend.".format(backend)
        )


def md(geom, backend="xtb", **kwargs):
    """
    Optimize geometry via molecular dyanmics using supplied forcefield
    and basis set.

    Parameters
    ----------
    backend : str
        Alias for backend selection (xtb).

    Returns
    -------
    result
        Object containing relevant outputs from the simulation.

    """

    # Select backend
    return _backend_selector(backend).run(geom, **kwargs)


class XTBWrapper(WrapperInterface):
    """
    Wrapper for xtb functionality.

    Implements :class:`~isicle.interfaces.WrapperInterface` to ensure required methods are exposed.

    Attributes
    ----------
    temp_dir : str
        Path to temporary directory used for simulation.
    geom : :obj:`isicle.geometry.Geometry`
        Internal molecule representation.
    result : dict
        Dictionary containing simulation results.

    """

    _defaults = ["geom", "result", "temp_dir"]
    _default_value = None

    def __init__(self):
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))
        self.temp_dir = isicle.utils.mkdtemp()

    def set_geometry(self, geom):
        """
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        """

        # Assign geometry
        self.geom = geom

        # Save geometry
        self.save_geometry()

    def save_geometry(self, fmt="xyz"):
        """
        Save internal :obj:`~isicle.geometry.Geometry` representation to file.

        Parameters
        ----------
        fmt : str
            Filetype used by xtb. Must be "xyz", "smi", ".inchi", ".mol", ".xyz",
            ".pdb", ".pkl".

        """
        # Path operations
        self._fmt = fmt.lower()
        geomfile = os.path.join(
            self.temp_dir, "{}.{}".format(self.geom.basename, self._fmt.lower())
        )

        # All other formats
        isicle.save(geomfile, self.geom)
        self.geom.path = geomfile

    def _configure_xtb(self, forcefield="gfn2", optlevel="normal", solvation=None):
        """
        Set command line for xtb simulations.

        Parameters
        ----------
        forcefield : str
            GFN forcefield for the optimization
            Default: gff
            Supported forcefields: gfn2, gfn1, gff
        optlevel : str
            Optimization convergence level
            Default : normal
            Supported : crude, sloppy, loose, lax, normal, tight, vtight extreme

        """

        # Add base command
        s = "xtb "

        # Add geometry
        s += "{}.{}".format(self.geom.basename, self._fmt.lower())

        # Add optimize tag
        s += " --opt " + optlevel + " "

        # Add forcefield
        s += "--" + forcefield + " "

        # Add charge
        s += "--chrg " + str(self.geom.formal_charge) + " "

        # Add optional implicit solvation
        if solvation is not None:
            s += "--alpb " + solvation + " "

        # Add output
        s += "&>" + " "

        s += "{}.{}".format(self.geom.basename, "out")
        return s

    def _configure_crest(
        self,
        ewin=6,
        optlevel="Normal",
        forcefield="gfn2",
        protonate=False,
        deprotonate=False,
        tautomerize=False,
        ion=None,
        dryrun=False,
        processes=1,
        solvation=None,
        ignore_topology=False,
    ):
        """
        Set command line for crest simulations.

        Parameters
        ----------
        ewin : int
            Energy window (kcal/mol) for conformer, (de)protomer, or tautomer search.
            Default : 6
        optlevel : str
            Optimization convergence level
            Default : normal
            Supported : crude, sloppy, loose, lax, normal, tight, vtight extreme
        forcefield : str
            GFN forcefield for the optimization
            Default: gff
            Supported forcefields: gfn2, gfn1, gff
        protonate : bool
            Signal to initiate protomer search. Suggested ewin = 30.
            Default : False
        deprotonate : bool
            Signal to initiate deprotonated conformers. Suggesting ewin = 30.
            Default : False
        tautomer : bool
            Signal to initiate tautomer search.
            Default : False
        ion : str
            Keyword to couple with protonate to ionize molecule with an ion other than a proton.
            See :obj:`~isicle.adduct.parse_ion` for list of ion options.

        """

        # Start base command
        s = "crest "

        # Add geometry
        s += str(
            os.path.join(
                self.temp_dir, "{}.{}".format(self.geom.basename, self._fmt.lower())
            )
        )

        s += " "
        # Add optional tag
        if protonate:
            s += "-protonate "
        elif deprotonate:
            s += "-deprotonate "
        elif tautomerize:
            s += "-tautomerize "

        if ion is not None:
            s += "-swel " + ion + " "

        # Add charge
        s += "-chrg " + str(self.geom.formal_charge) + " "

        # Add dryrun option
        if dryrun:
            s += "--dryrun "

        # Add energy window
        s += "--ewin " + str(ewin) + " "

        # Add optlevel
        s += "--optlevel " + optlevel + " "

        # Add forcefield
        s += "-" + forcefield + " "

        # Add optional solvation
        if solvation is not None:
            s += "--alpb " + solvation + " "

        # Number of processes
        s += "-T " + str(processes) + " "

        if ignore_topology:
            s += "--noreftopo "

        # Add output
        s += "&>" + " "

        s += os.path.join(self.temp_dir, "{}.{}".format(self.geom.basename, "out"))

        return s

    def configure(
        self,
        task="optimize",
        forcefield="gfn2",
        ewin=6,
        ion=None,
        optlevel="Normal",
        dryrun=False,
        processes=1,
        solvation=None,
        ignore_topology=False,
    ):
        """
        Generate command line

        Parameters
        ----------
        tasks : str
            Set task to "optimize", "conformer", "protonate", "deprotonate", or "tautomerize".
            Default : "optimize"
        forcefield : str
            GFN forcefield for the optimization
            Default: gff
            Supported forcefields: gfn2, gfn1, gff
        ewin : int
            Energy window (kcal/mol) for conformer(set to 6), (de)protomer(set to 30), or tautomer(set to 30) search.
            Default : 6
        ion : str
            Ion for protomer calculation.
        optlevel : str or list of str
            Set optimization level. Supply globally or per task.
        ion : str
            Keyword to couple with protonate to ionize molecule with an ion other than a proton.
            See :obj:`~isicle.adduct.parse_ion` for list of ion options.

        """

        if isinstance(task, list):
            raise TypeError("Initiate one xtb or crest job at a time.")
        if isinstance(forcefield, list):
            raise TypeError("Initiate one forcefield at a time.")
        if isinstance(optlevel, list):
            raise TypeError("Initiate one opt level at a time.")

        if task == "optimize":
            config = self._configure_xtb(optlevel=optlevel, forcefield=forcefield)

        else:
            if task == "conformer":
                p, d, t, i = False, False, False, None

            elif task == "protonate":
                p, d, t, i = True, False, False, ion

            elif task == "deprotonate":
                p, d, t, i = False, True, False, ion

            elif task == "tautomerize":
                p, d, t, i = False, False, True, ion

            if p is not None:
                config = self._configure_crest(
                    ewin=ewin,
                    optlevel=optlevel,
                    forcefield=forcefield,
                    protonate=p,
                    deprotonate=d,
                    tautomerize=t,
                    ion=i,
                    dryrun=dryrun,
                    processes=processes,
                    solvation=solvation,
                    ignore_topology=ignore_topology,
                )
            else:
                raise ValueError(
                    "Task not assigned properly, please choose optimize, conformer, protonate, deprotonate, or tautomerize"
                )

        self._task = task

        self._config = config

    def submit(self):
        """
        Run xtb or crest simulation according to configured inputs.
        """
        cwd = os.getcwd()
        os.chdir(self.temp_dir)
        subprocess.call(self._config, shell=True)
        os.chdir(cwd)

    def finish(self):
        """
        Parse results, save xtb output files, and clean temporary directory
        """

        # Get list of outputs
        outfiles = glob.glob(os.path.join(self.temp_dir, "*"))

        # Filter directories
        outfiles = [x for x in outfiles if not os.path.isdir(x)]

        # Result container
        result = {}

        # Split out geometry files
        geomfiles = sorted([x for x in outfiles if x.endswith(".xyz")])
        outfiles = sorted([x for x in outfiles if not x.endswith(".xyz")])

        # # Atom count
        # n_atoms = xtbw.geom.get_natoms()

        # Charge lookup
        charge_lookup = defaultdict(lambda: 0, {"protonated": 1, "deprotonated": -1})

        # Enumerate geometry files
        for geomfile in geomfiles:
            # Extract unique basename
            # Replace ancillary "crest_"
            basename = os.path.splitext(
                os.path.basename(geomfile).replace("crest_", "")
            )[0]

            # Only process of-interest structures
            if basename in [
                "struc",
                "best",
                "xtbopt",
                "protonated",
                "deprotonated",
                "tautomers",
                "conformers",
                "rotamers"
            ]:
                # Open file
                with open(geomfile, "r") as f:
                    contents = f.readlines()

                # Get file-specific atom count
                n_atoms = int(contents[0].strip())

                # Split into individual xyz
                contents = [
                    "".join(contents[i : i + n_atoms + 2])
                    for i in range(0, len(contents), n_atoms + 2)
                ]

                # Iterate XYZ content
                geoms = []
                for xyzblock in contents:
                    # Construct mol from XYZ
                    raw_mol = Chem.MolFromXYZBlock(xyzblock)
                    mol = Chem.Mol(raw_mol)
                    # rdDetermineBonds.DetermineBonds(
                    #     mol, charge=self.geom.formal_charge + charge_lookup[basename]
                    # )

                    # Initialize Geometry instance
                    geom = isicle.geometry.Geometry(
                        mol=mol, basename=self.geom.basename
                    )

                    # Update coordinates of starting geometry
                    geom = self.geom.update_coordinates(geom)

                    # Append to list of geometries
                    geoms.append(geom)

                # Add to result container
                if len(geoms) > 1:
                    result[basename] = geoms
                else:
                    result[basename] = geoms[0]

        # Enumerate output files
        for outfile in outfiles:
            # Split name and extension
            if "." in outfile:
                basename, ext = os.path.basename(outfile).rsplit(".", 1)

            # No extension
            else:
                basename = os.path.basename(outfile)
                ext = None

            # Key management
            if ext in [None, "tmp", "0"]:
                key = basename
            else:
                key = ext

            # Read output content
            with open(outfile, "rb") as f:
                contents = f.read()

            # Attempt utf-8 decode
            try:
                result[key] = contents.decode("utf-8")
            except UnicodeDecodeError:
                result[key] = contents

        # Renaming key
        rename = {
            "struc": "input",
            "xtbopt": "final",
            "original": "coord_original",
            "protonated": "cations",
            "deprotonated": "anions",
        }

        # Rename matching keys
        for key in result.keys() & rename.keys():
            result[rename[key]] = result.pop(key)

        # Assign to attribute
        self.result = result
        return self.result

    def parse(self):
        """
        Parse xTB simulation results.

        Returns
        -------
        dict
            Dictionary containing parsed outputs from the simulation.

        """

        if self.result is None:
            raise RuntimeError("Must complete xTB simulation.")

        parser = isicle.parse.XTBParser(data=self.result)

        return parser.parse()

    def run(self, geom, **kwargs):
        """
        Optimize geometry via density functional theory using supplied functional
        and basis set.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        **kwargs
            Keyword arguments to configure the simulation.
            See :meth:`~isicle.md.XTBWrapper.configure`.

        Returns
        -------
        :obj:`~isicle.md.XTBWrapper`
            Wrapper object containing relevant outputs from the simulation.

        """

        # New instance
        self = XTBWrapper()

        # Set geometry
        self.set_geometry(geom)

        # Configure
        self.configure(**kwargs)

        # Run QM simulation
        self.submit()

        # Finish/clean up
        self.finish()

        return self


class RDKitWrapper(Geometry, WrapperInterface):
    """
    Wrapper for RDKit functionality.

    Implements :class:`~isicle.interfaces.WrapperInterface` to ensure required methods are exposed.

    Attributes
    ----------
    geom : :obj:`isicle.geometry.Geometry`
        Internal molecule representation.
    method : str
        Method of RDKit conformer generation specified.
    numConfs : int
        The number of conformers to generate.
    result : dict
        Dictionary containing simulation results.
    """

    _defaults = ["geom", "method", "numConfs", "result"]
    _default_value = None

    def __init__(self):
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))

    def set_geometry(self, geom):
        """
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        """

        # Assign geometry
        self.geom = geom

    def configure(self, method: str = "etkdgv3", numConfs: int = 10, **kwargs):
        """
        Set conformer generation parameters.
        Parameters
        ----------
        method: str
            `distance` for distance geometry method (default)
            `etkdg` for ETKDG method
            `etkdgv2` for ETKDG method
            `etkdgv3` for ETKDG method
            `sretkdgv3` for ETKDG method
        numConfs: int
            the number of conformers to generate
        **kwargs:
            Keyword arguments to configure the simulation
            See :meth:`~isicle.md.RDKitWrapper._configure_distance_geometry`
        """
        lookup = {
            "distance": self._configure_distance_geometry,
            "etdg": self._configure_etdg,
            "etkdg": self._configure_etkdg1,
            "etkdgv1": self._configure_etkdg1,
            "etkdgv2": self._configure_etkdg2,
            "etkdgv3": self._configure_etkdg3,
            "sretkdgv3": self._configure_etkdg3_variant,
        }
        method = str(method)
        try:
            lookup[method.lower()](**kwargs)
        except KeyError:
            raise
        self.method = method.lower()
        try:
            numConfs = int(numConfs)
        except ValueError:
            raise
        self.numConfs = numConfs

    def _configure_distance_geometry(
        self,
        pruneRmsThresh: float = -1.0,
        forceTol: float = 0.001,
        randomSeed: int = -1,
    ):
        """
        Set parameters for distance geometry based conformer generation.
        Parameters
        ----------
        pruneRmsThresh: float
            Greedy pruning mainting conformers are <float> apart based upon RMSD of heavy atoms.
            Default: no pruning, -1.0
        forceTol: float
            Tolerance to be used in force-field minimizations
        randomSeed:  int
            provide a seed for the random number generator
            `-1` causes RNG to not be seeded
        """
        self.pruneRmsThresh = pruneRmsThresh
        self.forceTol = forceTol
        self.randomSeed = randomSeed

    def _configure_etdg(self):
        """
        Set parameters for ETDG conformer generation.
        """
        self.params = rdDistGeom.ETDG()

    def _configure_etkdg1(self):
        """
        Set parameters for ETKDG conformer generation, based on work by Riniker and Landrum.
        Version 1: RDKit default
        """
        self.params = rdDistGeom.ETKDG()

    def _configure_etkdg2(self):
        """
        Set parameters for ETKDG conformer generation, based on work by Riniker and Landrum.
        Version 2: (default) 2016 release, updated torsion angle potentials
        """
        self.params = rdDistGeom.ETKDGv2()

    def _configure_etkdg3(self, variant=True):
        """
        Set parameters for ETKDG conformer generation, based on work by Riniker and Landrum.
        Version 3: Updated sampling small rings AND macrocycles
        """
        self.params = rdDistGeom.ETKDGv3()
        if variant is True:
            self.params.useSmallRingTorsions = True

    def _configure_etkdg3_variant(self):
        """
        Set parameters for ETKDG conformer generation, based on work by Riniker and Landrum.
        Version 3 variant: Updated sampling for small rings, NOT macrocycles
        """
        self.params = rdDistGeom.srETKDGv3()

    def submit(self):
        """
        Execute conformer generation with user-specifed method, parameters.
        """
        if self.method == "distance":
            rdDistGeom.EmbedMultipleConfs(
                self.geom.mol,
                numConfs=self.numConfs,
                randomSeed=self.randomSeed,
                pruneRmsThresh=self.pruneRmsThresh,
                forceTol=self.forceTol,
            )
        elif "etkdg" in self.method or "etdg" in self.method:
            rdDistGeom.EmbedMultipleConfs(
                self.geom.mol, numConfs=self.numConfs, params=self.params
            )
        else:
            raise ValueError(
                "Failure to run RDKit MD, method and/or variant not recognized"
            )
    
    def parse(self):
        print("No need to parse RDKit result. Simply access `result` attribute.")

    def finish(self):
        """
        Parse RDKit conformers generated.
        """
        conf_count = self.geom.mol.GetNumConformers()
        conformers = [
            isicle.load(Chem.Mol(self.geom.mol, confId=i)) for i in range(conf_count)
        ]

        for conf, label in zip(conformers, range(conf_count)):
            conf.__dict__.update(conformerID=label, basename=self.geom.basename)

        self.result = isicle.conformers.ConformationalEnsemble(conformers)
        return self.result

    def run(self, geom, **kwargs):
        """
        Generate conformers using RKDit and supplied parameters.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        **kwargs
            Keyword arguments to configure the simulation.
            See :meth:`~isicle.md.RDKitWrapper.configure`.

        Returns
        -------
        :obj:`~isicle.md.RDKitWrapper`
            Wrapper object containing relevant outputs from the simulation.

        """

        # New instance
        self = RDKitWrapper()

        # Set geometry
        self.set_geometry(geom)

        # Configure
        self.configure(**kwargs)

        # Run QM simulation
        self.submit()

        # Finish/clean up
        self.finish()

        return self
