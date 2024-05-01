from collections import defaultdict
import glob
import os
import subprocess

from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds, rdDistGeom

import isicle
from isicle.geometry import Geometry
from isicle.interfaces import WrapperInterface
from isicle.parse import TINKERParser, XTBParser

"""
Files resulting from an xtb job always run in the same directory that the command is
issued in, no matter where the input is. Can direct the .log file, but no other files.
"""


def _program_selector(program):
    """
    Selects a supported molecular dynamics program for associated simulation.
    Currently only NWChem has been implemented.

    Parameters
    ----------
    program : str
        Alias for program selection (xtb).

    Returns
    -------
    program
        Wrapped functionality of the selected program. Must implement
        :class:`~isicle.interfaces.WrapperInterface`.

    """

    program_map = {"xtb": XTBWrapper, "rdkit": RDKitWrapper, "tinker": TINKERWrapper}

    if program.lower() in program_map.keys():
        return program_map[program.lower()]()
    else:
        raise ValueError(
            "{} not a supported molecular dynamics program.".format(program)
        )


def md(geom, program="xtb", **kwargs):
    """
    Optimize geometry via molecular dyanmics using supplied forcefield
    and basis set.

    Parameters
    ----------
    program : str
        Alias for program selection (xtb).

    Returns
    -------
    result
        Object containing relevant outputs from the simulation.

    """

    # Select program
    return _program_selector(program).run(geom, **kwargs)


class XTBWrapper(WrapperInterface):
    """
    Wrapper for xtb functionality.

    Implements :class:`~isicle.interfaces.WrapperInterface` to ensure required methods are exposed.

    Attributes
    ----------
    temp_dir : str
        Path to temporary directory used for simulation.
    task_map : dict
        Alias mapper for supported molecular dynamic presets. Includes
        "optimize", "crest", "nmr", "protonate", "deprotonate", and "tautomer".
    geom : :obj:`isicle.geometry.Geometry`
        Internal molecule representation.
    fmt : str
        File extension indicator.
    job_list : str
        List of commands for simulation.

    """

    _defaults = ["geom"]
    _default_value = None

    def __init__(self, **kwargs):
        self.temp_dir = isicle.utils.mkdtemp()
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))
        self.__dict__.update(**kwargs)

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
        self.fmt = fmt.lower()
        geomfile = os.path.join(
            self.temp_dir, "{}.{}".format(self.geom.basename, self.fmt.lower())
        )

        # All other formats
        isicle.io.save(geomfile, self.geom)
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
        s += "{}.{}".format(self.geom.basename, self.fmt.lower())

        # Add optimize tag
        s += " --opt " + optlevel + " "

        # Add forcefield
        s += "--" + forcefield + " "

        # Add charge
        s += "--chrg " + str(self.geom.get_charge()) + " "

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
                self.temp_dir, "{}.{}".format(self.geom.basename, self.fmt.lower())
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
        s += "-chrg " + str(self.geom.get_charge()) + " "

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

        self.task = task

        self.config = config

    def submit(self):
        """
        Run xtb or crest simulation according to configured inputs.
        """
        cwd = os.getcwd()
        os.chdir(self.temp_dir)
        subprocess.call(self.config, shell=True)
        os.chdir(cwd)

    def finish(self):
        """
        Parse results, save xtb output files, and clean temporary directory
        """

        # Get list of outputs
        outfiles = glob.glob(os.path.join(self.temp_dir, "*"))

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
                    rdDetermineBonds.DetermineBonds(
                        mol, charge=self.geom.get_charge() + charge_lookup[basename]
                    )

                    # Initialize Geometry instance
                    geom = isicle.geometry.Geometry(
                        mol=mol, basename=self.geom.basename
                    )

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

        return result

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
        result = self.finish()

        return result


class RDKitWrapper(Geometry, WrapperInterface):
    """
    Wrapper for RDKit functionality.

    Implements :class:`~isicle.interfaces.WrapperInterface` to ensure required methods are exposed.

    Attributes
    ----------
    geom : :obj:`isicle.geometry.Geometry`
        Internal molecule representation.
    method: str
        Method of RDKit conformer generation specified.
    numConfs: int
        The number of conformers to generate.
    """

    _defaults = ["geom", "method", "numConfs"]
    _default_value = None

    def __init__(self, **kwargs):
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))
        self.__dict__.update(**kwargs)

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

        return isicle.conformers.ConformationalEnsemble(conformers)

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


class TINKERWrapper(Geometry, WrapperInterface):
    """
    Wrapper for TINKER functionality.

    Implements :class:`~isicle.interfaces.WrapperInterface` to ensure required methods are exposed.

    Attributes
    ----------
    temp_dir : str
        Path to temporary directory used for simulation.
    task_map : dict
        Alias mapper for supported molecular dynamic presets. Includes
        "optimize", "conformer", "nmr", "protonate", "deprotonate", and "tautomer".
    geom : :obj:`isicle.geometry.Geometry`
        Internal molecule representation.
    fmt : str
        File extension indicator.
    job_list : str
        List of commands for simulation.

    """

    _defaults = ["geom"]
    _default_value = None

    def __init__(self, **kwargs):
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))
        self.__dict__.update(**kwargs)

    def _convert_to_tinkerxyz(self):
        """
        Convert mol to TINKER XYZ format using code from DP5, Goodman Lab.

        Parameters
        ----------
        geom : :obj: `~isicle.geometry.Geometry`
            Molecule representation.

        """

        # getting MMFF values for large atom types
        def getMMFF_large_atom_type(mmff_props, atom, m):
            small_to_large_list = [
                [[]],
                [[1]],
                [[3, "C"], [2, "C=C"]],
                [
                    [4, "C=O"],
                    [5, "C=N"],
                    [6, "NC(N)=N"],
                    [7, "CC=O"],
                    [8, "NC=O"],
                    [10, "NC(=O)N"],
                    [11, "OC=O"],
                    [12, "NC(=O)O"],
                    [13, "NC(=O)O"],
                    [14, "OC(=O)O"],
                    [15, "SC=O"],
                    [16, "NC=S"],
                    [17, "C=S(O)O"],
                    [18, "C=S=O"],
                    [19, "SC=S"],
                    [20, "C=P"],
                ],
                [[21, "C#[C,N]"], [22, "[C,N,O]=C=[C,N,O]"]],
                [[23, "C[H]"], [24, "[Si][H]"]],
                [
                    [41, "O"],
                    [25, "OC"],
                    [26, "OC=O"],
                    [27, "OC=C"],
                    [27, "Occ"],
                    [28, "OC=N"],
                    [29, "OC=S"],
                    [31, "ON=O"],
                    [30, "O[N+]([O-])=O"],
                    [36, "OS"],
                    [34, "OSO"],
                    [35, "OS=O"],
                    [33, "OS(O)=O"],
                    [32, "OS(O)(=O)=O"],
                    [40, "OP"],
                    [39, "OPO"],
                    [38, "OP(=O)O"],
                    [37, "OP(=O)(=O)O"],
                ],
                [
                    [42, "C=O"],
                    [44, "CC=O"],
                    [43, "NC=O"],
                    [45, "OC=O"],
                    [46, "O=N"],
                    [47, "S=O"],
                    [48, "[C,N]=S=O"],
                ],
                [[49]],
                [[50, "C=N"], [51, "N=N"]],
                [[52, "NC=O"], [53, "NC=S"], [54, "NN=C"], [55, "NN=N"]],
                [[56]],
                [[57]],
                [[58]],
                [[59]],
                [[60]],
                [[61]],
                [[62, "S=O"], [63, "S=N"]],
                [
                    [64, "O=S=O"],
                    [70, "OSN"],
                    [65, "N-S(=O)=O"],
                    [66, "OS(O)O"],
                    [67, "C"],
                    [68, "OS(O)(O)O"],
                    [69, "CS(O)(O)C"],
                ],
                [[71]],
                [[72]],
                [[74, "[H]O"], [73, "[H]OC"], [75, "[H][O-]"]],
                [[76]],
                [
                    [82, "[H]N"],
                    [77, "[H]N(C)C"],
                    [78, "[H]N([H])[H]"],
                    [79, "[H]n1cccc1"],
                    [80, "[H]NO"],
                    [81, "[H][N-]"],
                ],
                [[83, "[H]OC=O"], [84, "[H]OP"]],
                [[89, "P"], [88, "PO"], [87, "OPO"], [86, "OP(O)O"], [89, "OP(O)(O)O"]],
                [[90]],
                [[91, "[H]N=N"], [92, "[H]N=C"]],
                [
                    [102, "[H]N"],
                    [93, "[H]NC=O"],
                    [94, "[H]NC=S"],
                    [95, "[H]NC=C"],
                    [96, "[H]NC=N"],
                    [97, "[H]NN=C"],
                    [98, "[H]NN=N"],
                    [99, "[H]NS=O"],
                    [100, "[H]NP=O"],
                    [101, "[H]N#[C,N]"],
                ],
                [[103, "[H]OC=C"], [104, "[C]OC=N"]],
                [[105]],
                [[106]],
                [
                    [107, "[O-]C=O"],
                    [108, "NO"],
                    [109, "ON=O"],
                    [110, "O[N+]([O-])=O"],
                    [111, "[O-][N+]([O-])=O"],
                    [112, "OS"],
                    [113, "OS=O"],
                    [114, "OS(=O)=O"],
                    [115, "OS(=O)(=O)O"],
                    [116, "O=[S-]S"],
                    [117, "OP"],
                    [118, "OPO"],
                    [119, "OP(=O)O"],
                    [120, "OP(=O)(=O)"],
                    [121, "OCl(=O)(=O)[O-]"],
                ],
                [[122]],
                [[123]],
                [[124, "[O-]"], [125, "[O-]C=[C,N]"]],
                [
                    [126, "[H][N+][H,C][H,C][H,C]"],
                    [127, "C1=[NH+]C=CN1"],
                    [128, "C1=C[NH+]=CC=C1"],
                    [129, "CC(N)=[NH2+]"],
                    [130, "C=[NH2+]"],
                    [131, "NC(N)=[NH2+]"],
                    [132, "[H]N([H])([H])([H])[H]"],
                ],
                [[133]],
                [[134]],
                [[135]],
                [[136, "NC=C"], [137, "NC=N"], [138, "NC=P"], [139, "NC#C"]],
                [[140, "[O-]C=O"], [141, "[S-]C=S"]],
                [[142]],
                [
                    [143, "NS(=O)O"],
                    [144, "NS(=O)(=O)O"],
                    [145, "NP(=O)O"],
                    [146, "NP(=O)(=O)O"],
                    [147, "NC#N"],
                ],
                [[148]],
                [[149, "ON=O"], [150, "O[N+][O-]=O"]],
                [[151]],
                [[152]],
                [[153]],
                [[154]],
                [[155]],
                [[156]],
                [[157]],
                [[158]],
                [[159, "[N+]=C"], [160, "[N+]=N"]],
                [[161]],
                [[162]],
                [[163, "NC(N)=[NH2+]"], [164, "[N+]=CN"]],
                [[165]],
                [[166]],
                [[167]],
                [[168]],
                [[169]],
                [[170]],
                [[171]],
                [[172]],
                [[173]],
                [[174]],
                [[175]],
                [[176]],
                [[177]],
                [[178, "[H]S"], [179, "[H]S=N"], [180, "[H]P"]],
                [[181, "SP"], [183, "[S-]"], [182, "[S-]C=S"], [184, "[S-]S(=O)"]],
                [[185, "[O-]S=O"], [186, "[O-]S=S"]],
                [[187]],
                [[188]],
                [[189]],
                [[190]],
                [[191]],
                [[192]],
                [[193]],
                [
                    [194, "N[N+]1=CNC=C1"],
                    [195, "[H][N+]([H])([H])([H])[H]"],
                    [196, "[H][N+]([H])([H])([H])[H]"],
                    [197, "[H][N+]([H])([H])([H])[H]"],
                ],
                [
                    [198, "[H][N+]([H])([H])([H])[H]"],
                    [199, "[H][N+]([H])([H])([H])[H]"],
                    [200, "[H][N+]([H])([H])([H])[H]"],
                ],
                [[-1]],
                [[-1]],
                [[-1]],
                [[-1]],
                [[201]],
                [[202]],
                [[203]],
                [[204]],
                [[205]],
                [[206]],
                [[207]],
                [[208]],
                [[209, "[Zn]"], [210, "[Zn++]"]],
                [[211]],
                [[212]],
                [[213]],
                [[214]],
            ]

            MMFF_small_atom_type = mmff_props.GetMMFFAtomType(atom.GetIdx())
            MMFF_large_atom_type = small_to_large_list[MMFF_small_atom_type][0][0]
            if len(small_to_large_list[MMFF_small_atom_type]) > 1:
                for atom_info in small_to_large_list[MMFF_small_atom_type]:
                    substructure = Chem.MolFromSmarts(atom_info[1])
                    for substructure_match in m.GetSubstructMatches(substructure):
                        if substructure_match.count(atom.GetIdx()) > 0:
                            MMFF_large_atom_type = atom_info[0]
            return MMFF_large_atom_type

        mol = self.geom.mol
        conf = mol.GetConformer()
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol)

        xyz = ""

        # Set header line with number of atoms and basename
        xyz += "{:>6}  {}\n".format(mol.GetNumAtoms(), self.geom.basename)

        for atom in mol.GetAtoms():
            bond_list = []
            attached_atoms = ""

            for connection in atom.GetNeighbors():
                bond_list.append(connection.GetIdx() + 1)
            bond_list.sort()

            for connection in bond_list:
                attached_atoms += "{:>5}".format(str(connection)) + " "

            xyz += "{:>6} {:>2} {:13.6f} {:11.6f} {:11.6f} {:>5} {}\n".format(
                atom.GetIdx() + 1,
                atom.GetSymbol(),
                list(conf.GetAtomPosition(atom.GetIdx()))[0],
                list(conf.GetAtomPosition(atom.GetIdx()))[1],
                list(conf.GetAtomPosition(atom.GetIdx()))[2],
                getMMFF_large_atom_type(mmff_props, atom, mol),
                attached_atoms,
            )

        return xyz

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

        self.tinkerxyz = self._convert_to_tinkerxyz()

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
        # Path operationspyth
        self.temp_dir = isicle.utils.mkdtemp()
        self.fmt = fmt.lower()
        geomfile = os.path.join(
            self.temp_dir, "{}.{}".format(self.geom.basename, self.fmt.lower())
        )

        with open(geomfile, "w+") as f:
            f.write(self.tinkerxyz)
        f.close()

        self.geom.path = geomfile

    def configure(self, task="scan", tinker_path="~/tinker"):
        config = tinker_path + "/bin/" + task + " " + self.geom.path + " "
        config += tinker_path + "/params/mmff.prm 0 10 20 0.00001 "
        config += "| tee ./" + self.geom.basename + ".tout"

        self.config = config

    def submit(self):
        owd = os.getcwd()
        os.chdir(self.temp_dir)
        job = self.config
        subprocess.call(job, shell=True)
        os.chdir(owd)

    def finish(self):
        parser = TINKERParser()

        parser.load(os.path.join(self.temp_dir, self.geom.basename + ".tout"))
        self.output = parser.load(
            os.path.join(self.temp_dir, self.geom.basename + ".tout")
        )

        result = parser.parse()

        self.__dict__.update(result)

        conformerID = 1
        for i in self.geom:
            i.add___dict__({k: v for k, v in result.items() if k != "geom"})
            i.__dict__.update(basename=self.geom.basename)
            i.__dict__.update(conformerID=conformerID)
            conformerID += 1
            return self

        return self

    def run(self, geom, **kwargs):
        """
        Take geometry and conduct specified tasks using TINKER

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        **kwargs
            Keyword arguments to configure the simulation.
            See :meth:`~isicle.md.XTBWrapper.configure`.

        Returns
        -------
        :obj:`~isicle.md.TINKERWrapper`
            Wrapper object containing relevant outputs from the simulation.

        """

        # New instance
        self = TINKERWrapper()

        # Set geometry
        self.set_geometry(geom)

        # Configure
        self.configure(**kwargs)

        # Run QM simulation
        self.submit()

        # Finish/clean up
        self.finish()

        return self
