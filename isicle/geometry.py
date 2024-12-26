import numpy as np
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import AllChem

import isicle
from isicle.interfaces import GeometryInterface


class Geometry(GeometryInterface):
    """
    Molecule information, specialized for 3D representations. It is not
    recommended to manipulate or retrieve attributes of this class without
    using class functions.
    """

    _defaults = ["mol", "charge", "basename"]
    _default_value = None

    def __init__(self, **kwargs):
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))
        self.__dict__.update(kwargs)

        # if self.charge is None:
        #     self.charge = self.get_charge()

    def view(self):
        return self.to_mol()

    def get_basename(self):
        """
        Returns a copy of this object's basename (original filename).

        Returns
        -------
        str
            Basename of original filepath.

        """

        return self.basename

    def add___dict__(self, d, override=False):
        """
        Accepts a dictionary of values and adds any non-conflicting
        information to the attribute dictionary.

        Parameters
        ----------
        d : dict
            Attribute dictionary.
        override : bool
            Singal whether to override existing attributes.

        """

        if not override:
            # Remove keys that are already present in attribute dictionary
            [d.pop(k, None) for k in self.__dict__.keys()]

        # Add to attributes
        self.__dict__.update(d)

    def get___dict__(self):
        """
        Return a copy of this object's attributes as a dictionary.

        Returns
        -------
        dict
            Instance attributes.

        """

        return self.__dict__.copy()

    def __copy__(self, **kwargs):
        """
        Return hard copy of this class instance.

        Parameters
        ----------
        kwargs
            Keyword arguments to update copy.

        Returns
        -------
        :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        """

        # Copy attributes
        d = self.__dict__.copy()

        # Safely copy mol if present
        if "mol" in d:
            d["mol"] = self.to_mol()

        # Build self instance from scratch
        instance = type(self)(**d)

        # Update from kwargs
        instance.__dict__.update(kwargs)

        return instance

    def _is_embedded(self, mol):
        """
        Check if molecule has embedded 3D coordinates.

        Returns
        -------
        bool
            Indication of 3D coordinate presence.

        """

        try:
            mol.GetConformer().Is3D()
            return True
        except:
            return False

    def addHs(self):
        """
        Add implicit hydrogens to molecule.

        Returns
        -------
         :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        """

        # Get copy of mol object
        mol = self.to_mol()

        # Add Hs with coordinates
        mol = Chem.AddHs(mol, addCoords=True)

        # Return new geometry instance
        return self.__copy__(mol=mol)

    def initial_optimize(self, embed=False, forcefield="UFF", ff_iter=200):
        """
        Initial molecule optimization by basic force field to establish rough
        3D coordinates. Further optimization (molecular dynamics, density
        functional theory) recommended.

        Parameters
        ----------
        embed: bool
            Attempt molecule embedding with eigenvalues of distance matrix.
            Failure results in seeding with random coordinates.
        forcefield: str
            Specify "UFF" for universal force field, "MMFF" or "MMFF94" for
            Merck molecular force field 94, or "MMFF94s" for the MMFF94 s variant.

        Returns
        -------
        :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        """

        def _forcefield_selector(forcefield, mw):
            """
            Checks if user specified forcefield is compatible with supplied mol
            """
            forcefield = forcefield.lower()
            if forcefield == "uff":
                if Chem.rdForceFieldHelpers.UFFHasAllMoleculeParams(mw) is True:
                    return AllChem.UFFOptimizeMolecule
                else:
                    raise ValueError("UFF is not available for all atoms in molecule.")
            elif forcefield in ["mmff", "mmff94", "mmff94s"]:
                if Chem.rdForceFieldHelpers.MMFFHasAllMoleculeParams(mw) is True:
                    return Chem.rdForceFieldHelpers.MMFFOptimizeMolecule
                else:
                    raise ValueError(
                        "specified MMFF is not available for all atoms in molecule."
                    )
            else:
                raise ValueError(
                    "RDKit only supports UFF, MMFF, MMFF94, MMFF94s as forcefields."
                )

        # Get rdkit mol
        mol = self.to_mol()

        # Embed molecule 3D coordinates
        if embed is True:
            # Attempt embedding
            res = AllChem.EmbedMolecule(mol)
            if res == -1:
                # Use random coordinates
                res = AllChem.EmbedMolecule(mol, useRandomCoords=True)
                if res == -1:
                    raise ValueError("Embedding failure.")

        # Optimize according to supplied forcefield
        if forcefield is not None:
            # Check if embedded
            if self._is_embedded(mol):
                # Forcefield selection
                if "mmff94s" in forcefield.lower():
                    _forcefield_selector(forcefield, mol)(
                        mol, mmffVariant=forcefield, maxIters=ff_iter
                    )
                else:
                    _forcefield_selector(forcefield, mol)(mol, maxIters=ff_iter)

            # Not embedded
            else:
                raise ValueError("Molecule must have embedded 3D coordinates.")

        # Return copy with updated mol
        return self.__copy__(mol=mol)

    def desalt(self, salts=None):
        """
        Desalts RDKit mol object using Chem.SaltRemover module.

        Parameters
        ----------
        salts : str (optional)
            Salt type to remove. Ex: 'Cl', 'Br', '[Na+]'. Default: None.

        Returns
        -------
        :obj:`~isicle.geometry.Geometry`
            Moleculerepresentation.

        """

        # TODO: should this raise an error instead?
        # If no salts given, skip desalting
        if salts is None:
            return self.__copy__()

        remover = SaltRemover(defnFormat="smiles", defnData=salts)
        # defnData="[Cl,Br,Na]" *sample definition of salts to be removed*
        # add iterator for salts listed in config?
        # set potential salts to be removed in a config file

        mol, deleted = remover.StripMolWithDeleted(self.to_mol())
        # using StripMolWithDeleted instead of StripMol
        # add functionality to track removed salts
        # atomno = res.GetNumAtoms
        # if relevant to future use, returns atom count post desalting

        return self.__copy__(mol=mol)

    def neutralize(self):
        """
        Neutralizes RDKit mol object using neutralization reactions.

        Returns
        -------
        :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        """

        def _initialize_neutralisation_reactions():
            patts = (
                # Imidazoles
                ("[n+;H]", "n"),
                # Amines
                ("[N+;!H0]", "N"),
                # Carboxylic acids and alcohols
                ("[$([O-]);!$([O-][#7])]", "O"),
                # Thiols
                ("[S-;X1]", "S"),
                # Sulfonamides
                ("[$([N-;X2]S(=O)=O)]", "N"),
                # Enamines
                ("[$([N-;X2][C,N]=C)]", "N"),
                # Tetrazoles
                ("[n-]", "[nH]"),
                # Sulfoxides
                ("[$([S-]=O)]", "S"),
                # Amides
                ("[$([N-]C=O)]", "N"),
            )
            return [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y)) for x, y in patts]

        reactions = _initialize_neutralisation_reactions()

        mol = self.to_mol()

        # TODO: `replaced` is unused
        replaced = False
        for i, (reactant, product) in enumerate(reactions):
            while mol.HasSubstructMatch(reactant):
                replaced = True
                rms = AllChem.ReplaceSubstructs(mol, reactant, product)
                mol = rms[0]

        return self.__copy__(mol=mol)

    def tautomerize(self, return_all=False):
        """
        Generate tautomers according to RDKit TautomerEnumerator() method.

        Parameters
        ----------
        return_all : boolean (optional)
            If true, return all tautomers generated. Otherwise, only return
            the most common. Default=False

        Returns
        -------
        :obj:`~isicle.geometry.Geometry` or list of :obj:`~isicle.geometry.Geometry`
            Tautomer(s) of starting structure.

        """

        # source: https://rdkit.blogspot.com/2020/01/trying-out-new-tautomer.html
        # Discuss noted double bond changes
        enumerator = rdMolStandardize.TautomerEnumerator()

        mol = self.to_mol()
        res = [mol]
        tauts = enumerator.Enumerate(mol)
        smis = [Chem.MolToSmiles(x) for x in tauts]
        s_smis = sorted((x, y) for x, y in zip(smis, tauts) if x != self.to_smiles())
        res += [y for x, y in s_smis]

        # Ensure res is a list of mol objects
        if return_all:
            new_geoms = []
            for r in res:
                geom = self.__copy__(mol=r)
                new_geoms.append(geom)
            return new_geoms

        return self.__copy__(mol=res[0])

    def kekulize(self):
        """
        Kekulizes the molecule if possible. Otherwise the molecule is not modified.
        This is recommended for aromatic molecules undergoing explicit ionization.
        Aromatic bonds are explicitly defined and aromatic flags are cleared.

        """

        mol = Chem.KekulizeIfPossible(self.to_mol(), clearAromaticFlags=True)
        return self.__copy__(mol=mol)

    def ionize(self, ion_path=None, ion_list=None, method="explicit", **kwargs):
        """
        Ionize geometry, using specified list of ions and method of ionization.

        Parameters
        ----------
        ion_path : str
            Filepath to text file containing ions with charge (eg. `H+`) to be
            considered. Either ion_path or ion_list must be specified.
        ion_list : list
            List of strings of adducts to be considered. Must be specifed in
            syntax `Atom+` or `Atom-`, eg. `H+`, `Na+`, `H-Na+`. Either ion_path
            or ion_list must be specified
        method : str
            Method of ionization to be used, 'explicit' or 'crest' is accepted
        ensembl : bool (optional)
            Returns instead a list of adduct geometries
        kwargs:
            see :meth: `~isicle.adducts.ExplicitIonizationWrapper.submit` or
                       `~isicle.adducts.CRESTIonizationWrapper.submit`
            for more options

        Returns
        -------
        :obj:`~isicle.adducts.IonizationWrapper`
            Ionization result.

        """

        iw = isicle.adducts.ionize(method).run(
            self.__copy__(), ion_path=ion_path, ion_list=ion_list, **kwargs
        )

        return iw

    def get_natoms(self):
        """
        Calculate total number of atoms.

        Returns
        -------
        int
            Number of atoms.

        """

        return Chem.Mol.GetNumAtoms(self.to_mol())

    def get_atom_indices(
        self, atoms=["C", "H"], lookup={"C": 6, "H": 1, "N": 7, "O": 8, "F": 9, "P": 15}
    ):
        """
        Extract indices of each atom from the internal geometry.

        Parameters
        ----------
        atoms : list of str
            Atom types of interest.
        lookup : dict
            Mapping between atom symbol and atomic number.

        Returns
        -------
        list of int
            Atom indices.

        """

        atoms = [lookup[x] for x in atoms]
        idx = []
        for a in self.mol.GetAtoms():
            if a.GetAtomicNum() in atoms:
                idx.append(a.GetIdx())

        return idx

    def get_total_partial_charge(self):
        """
        Sum the partial charge across all atoms.

        Returns
        -------
        float
            Total partial charge.

        """

        mol = self.to_mol()
        AllChem.ComputeGasteigerCharges(mol)
        contribs = [
            mol.GetAtomWithIdx(i).GetDoubleProp("_GasteigerCharge")
            for i in range(mol.GetNumAtoms())
        ]
        return np.nansum(contribs)

    def get_charge(self):
        """
        Get formal charge of the molecule.

        Returns
        -------
        int
            Formal charge.

        """

        return Chem.rdmolops.GetFormalCharge(self.to_mol())

    def dft(self, backend="NWChem", **kwargs):
        """
        Perform density functional theory calculations according to supplied task list
        and configuration parameters for specified backend.

        Parameters
        ----------
        backend : str
            Alias for backend selection (NWChem, ORCA).
        kwargs
            Keyword arguments supplied to selected program. See `run` method of relevant
            wrapper for configuration parameters, e.g. :meth:`~isicle.qm.NWChemWrapper.run`.

        Returns
        -------
        :obj:`~isicle.qm.{program}Wrapper`
            Wrapper object containing relevant outputs from the simulation.

        """

        return isicle.qm.dft(self, backend=backend, **kwargs)

    def md(self, program="xtb", **kwargs):
        """
        Optimize geometry or generate conformers or adducts using stated forcefield.

        Parameters
        ----------
        kwargs
            Keyword arguments supplied to selected program. See `run` method of relevant
            wrapper for configuration parameters, e.g. :meth:`~isicle.md.XTBWrapper.run`.

        Returns
        -------
        :obj:`~isicle.md.{program}Wrapper`
            Wrapper object containing relevant outputs from the simulation.

        """

        return isicle.md.md(self, program=program, **kwargs)

    def to_mol(self):
        """
        Returns :obj:`~rdkit.Chem.rdchem.Mol` instance for this Geometry.

        Returns
        -------
        :obj:`~rdkit.Chem.rdchem.Mol`
            RDKit Mol instance.

        """

        return self.mol.__copy__()

    def to_smiles(self):
        """
        Get SMILES for this structure.

        Returns
        -------
        str
            SMILES string.

        """

        return Chem.MolToSmiles(self.to_mol())

    def to_inchi(self):
        """
        Get InChI for this structure.

        Returns
        -------
        str
            InChI string.

        """

        return Chem.MolToInchi(self.to_mol())

    def to_smarts(self):
        """
        Get SMARTS for this structure.

        Returns
        -------
        str
            SMARTS string.

        """

        return Chem.MolToSmarts(self.to_mol())

    def to_xyzblock(self):
        """
        Get XYZ text for this structure.

        Returns
        -------
        str
            XYZ representation as string.

        """

        return Chem.MolToXYZBlock(self.to_mol())

    def to_pdbblock(self):
        """
        Get PDB text for this structure.

        Returns
        -------
        str
            PDB representation as string.

        """

        return Chem.MolToPDBBlock(self.to_mol())

    def to_molblock(self):
        """
        Get :obj:`~rdkit.Chem.rdchem.Mol` text for this structure.

        Returns
        -------
        str
            :obj:`~rdkit.Chem.rdchem.Mol` representation as string.

        """

        return Chem.MolToMolBlock(self.to_mol())
