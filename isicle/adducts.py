import functools
import operator
import os
import pickle
import re

import rdkit.Chem.rdchem as rdc
from rdkit import Chem

import isicle
from isicle.conformers import ConformationalEnsemble
from isicle.interfaces import WrapperInterface
from isicle.md import md
from isicle.utils import safelist


def _parse_file_ions(path):
    """
    Read adduct ions (eg. 'Na+', 'H+', 'Mg2+') from file.

    Parameters
    ----------
    path : str
        Path to ions text file (.txt, etc.), each ion should be a new line.

    Returns
    -------
    parsed_contents
        List of ions from given text file.
    """
    if path is None:
        raise TypeError("Path to file containing ions must be specified")
    if os.path.isfile(path) == True:
        contents = isicle.geometry._load_text(path)
        parsed_contents = [line.rstrip() for line in contents]
        return parsed_contents
    else:
        raise TypeError("Path to file is not valid")


def _parse_list_ions(ion_list):
    """
    Parse and categorize list of ions by ion type.

    Input
    -----
    ion_list: list
        List of all ions to be considered, eg. ['K+','H-','H-Na+']
        Strings in list must have notation 'IonCharge' eg. 'H+', not '+H'
        Complex ions should be arranged in order of desired operation.
        (i.e. 'H-K+' specifies deprotonatin first, potassium added second)

    Returns
    ----------
    anions : list
        List of anions, default None
    cations : list
        List of cations, default None
    complex : list
        List of complex ions, default None
    """
    anions = []
    cations = []
    complex = []
    ion_list = safelist(ion_list)
    for x in ion_list:
        if ("+" in x) and ("-" in x):
            complex.append(x)
        elif "+" in x:
            cations.append(x)
        elif "-" in x:
            anions.append(x)
        else:
            raise ValueError("Unrecognized ion specified: {}".format(x))
    return (cations, anions, complex)


def _parse_ions(ion_path=None, ion_list=None):
    """
    Parses supplied ion information to recognized format

    Params
    ------
    ion_path: str
        Filepath to text file with an ion on each new line
        For more layout information, see: adducts._parse_file_ions()
    ion_list: list of str
        List of ions in format 'ioncharge' (eg. 'H+' or 'Na+')
        For more information on complex ions, see: adducts._parse_list_ions()

    Returns
    -------
    Tuple of (cations, anions, complex ions)
    """
    if ion_path is not None:
        # Load ion file
        ion_list = _parse_file_ions(ion_path)
    if ion_list is not None:
        # Parse ion file
        return _parse_list_ions(ion_list)
    else:
        raise RuntimeError("No ions supplied to parse.")


def _check_atom_group(ion_atomic_num):
    """
    Checks periodic group atom belongs to

    Params
    ------
    ion_atomic_num: Int

    Returns
    -------
    Bool
    """
    return ion_atomic_num in [1, 3, 4, 11, 12, 19, 20, 37, 38, 55, 56, 87, 88]


def _filter_by_substructure_match(init_mol, unknown_valid_list):
    """
    Filters ions in supplied list that are a substructure match to the supplied mol

    Params
    ------
    init_mol: RDKit mol object
    unknown_valid_list: list of str
        List of ions (eg. [`H-`])
    """
    valid_list = []
    for ion in unknown_valid_list:
        parsed_ion = re.findall(".+?[0-9]?[+-]", ion)
        subset = [
            f"[{re.findall('(.+?)[0-9]?[-]', i)[0]}]" for i in parsed_ion if "-" in i
        ]
        subset_mol = [Chem.MolFromSmiles(i, sanitize=False) for i in subset]
        substruc_check = [i for i in subset_mol if init_mol.HasSubstructMatch(i)]
        if False in substruc_check:
            continue
        else:
            valid_list.append(ion)
    return valid_list


def get_nested_dict_values(d):
    for v in d.values():
        if isinstance(v, dict):
            yield from get_nested_dict_values(v)
        else:
            yield v


class AdductEnsemble(dict):
    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self):
        return repr(self.__dict__)

    def __len__(self):
        return len(self.__dict__)

    def __delitem__(self, key):
        del self.__dict__[key]

    def clear(self):
        return self.__dict__.clear()

    def copy(self):
        return self.__dict__.copy()

    def has_key(self, k):
        return k in self.__dict__

    def update(self, *args, **kwargs):
        return self.__dict__.update(*args, **kwargs)

    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()

    def get_structures(self):
        values = get_nested_dict_values(self.__dict__)
        return isicle.conformers.build_conformational_ensemble(values)


def build_adduct_ensemble(geometries):
    # Cast to conformational ensemble
    if not isinstance(geometries, isicle.conformers.ConformationalEnsemble):
        geometries = isicle.conformers.build_conformational_ensemble(geometries)

    geometries._check_attributes("ion")
    geometries._check_attributes("adductID")

    d = AdductEnsemble()
    for adduct in geometries:
        if adduct.ion not in d.keys():
            d[adduct.ion] = {}

        d[adduct.ion][adduct.adductID] = adduct

    return d


def ionize(ion_method):
    """
    Selects a supported ionization method.

    Parameters
    ----------
    ion_method : str
        Alias for ion method selection (e.g. explicit).

    Returns
    -------
    program
        Wrapped functionality of the selected program. Must implement
        :class:`~isicle.interfaces.AdductInterface`.
    """

    ion_method_map = {
        "explicit": ExplicitIonizationWrapper,
        "crest": CRESTIonizationWrapper,
    }

    if ion_method.lower() in ion_method_map.keys():
        return ion_method_map[ion_method.lower()]()
    else:
        raise ValueError("{} not a supported ionization method.".format(ion_method))


def write(IonizationWrapper, path, fmt):
    """
    Write mol objects to file.

    Parameters
    ----------
    path : str
        Directory to write output files
    fmt : str
        Format in which to save the RDKit mol object (eg. 'mol2', 'pdb')
    """

    if path is not None and fmt is not None:
        for adduct in IonizationWrapper.adducts:
            isicle.geometry.Geometry.save(
                adduct.mol,
                os.path.join(
                    path, "{}_{}".format(adduct.basename, adduct.ion, adduct.adductID)
                ),
                fmt,
            )
        return
    elif path is not None:
        raise (
            "path passed to isicle.adducts.write; fmt flag must also be passed; data not written"
        )
    elif fmt is not None:
        raise (
            "fmt is supplied to isicle.adducts.write; path flag must also be passed; data not saved"
        )
    else:
        return


class ExplicitIonizationWrapper(WrapperInterface):
    """ """

    _defaults = ("geom", "adducts")
    _default_value = None

    def __init__(self, **kwargs):
        """
        Initialize :obj:`~isicle.adducts.ExplicitIonizationWrapper` instance.
        """
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))
        self.__dict__.update(**kwargs)
        if self.adducts is None:
            self.adducts = {}

    def set_geometry(self, geom):
        """
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.

        Parameters
        ----------
        geometry : :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        """

        if self.__dict__.get("geom") is not None:
            pass
        elif self.__dict__.get("mol") is not None:
            self.geom = geom
        else:
            raise ValueError("Could not find self.mol or self.geom")

    def set_charge(self):
        """ """

        if self.geom.__dict__.get("charge") is None:
            self.geom.__dict__.update(charge=self.geom.get_formal_charge())

    def _set_ions(self, ion_path=None, ion_list=None):
        """ """

        cations, anions, complex = _parse_ions(ion_path=ion_path, ion_list=ion_list)
        # Cast to list safely
        self.adducts["cations"] = safelist(cations)
        self.adducts["anions"] = safelist(anions)
        self.adducts["complex"] = safelist(complex)

    def _check_valid(self):
        """
        Perform substructure search.

        """

        self.adducts["anions"] = _filter_by_substructure_match(
            self.geom.mol, self.adducts["anions"]
        )
        self.adducts["complex"] = _filter_by_substructure_match(
            self.geom.mol, self.adducts["complex"]
        )

    def configure(self, ion_path=None, ion_list=None):
        """ """
        self._set_ions(ion_path=ion_path, ion_list=ion_list)
        self._check_valid()

    def _modify_atom_dict(
        self, init_mol, ion_atomic_num, mode=None, include_Alkali_ne=False
    ):
        """
        Downselect atom sites that can accept ions

        Parameters
        ----------
        init_mol: RDKit mol object
        ion_atomic_num: Int
            Integer denoting the atomic number of the ion to be added or removed
        mode: str
            'positive' for positive mode
            'negative' for negative mode
        include_Alkaline_ne: Bool
            If False: removes alkali and alkaline metals from base atoms to be ionized (default)
            If True: Adducts are attempted to be formed at these atoms, in addition to other atoms.

        Returns
        -------
        all_atom_dict defined as {<atom.rdkit.obj>:['C',3,2,1]}
            -atom.GetSymbol() returns atomic symbol e.g. 'C'
            -atom.GetTotalValence() returns total (explicit+implicit) valence
            -atom.GetDegree() returns number of directly-bonded neighbours
             indep. of bond order, dep. of Hs set to explicit
            -atom.GetTotalNumHs returns total (explicit+implicit) Hs on atom

        """

        pt = Chem.GetPeriodicTable()
        all_atom_dict = {
            atom.GetIdx(): [
                atom.GetSymbol(),
                atom.GetTotalValence() - atom.GetDegree(),
                atom.GetTotalNumHs(includeNeighbors=True),
                atom,
            ]
            for atom in init_mol.GetAtoms()
        }

        if mode == "positive":
            atom_dict = {}
            for key, value in all_atom_dict.items():
                # Check if atom atomic number is the same as the specified ion
                atom_atomic_num = pt.GetAtomicNumber(value[0])
                if include_Alkali_ne == False:
                    # Remove ion sites that are alkali or alkaline metals
                    if _check_atom_group(atom_atomic_num):
                        continue
                # Remove atoms that cannot accept a bond & do not have H to bump
                elif value[1] == 0 and value[2] == 0:
                    continue
                atom_dict[key] = value
            return atom_dict

        elif mode == "negative":
            atom_dict = {}
            for key, value in all_atom_dict.items():
                # Check if ion atomic num is in list of neighbouring atomic nums
                nbrs = [nbr.GetAtomicNum() for nbr in value[3].GetNeighbors()]
                if ion_atomic_num in nbrs:
                    atom_dict[key] = value
                else:
                    continue

            return atom_dict

        elif mode == None:
            raise ValueError("Must supply an ionization mode to _modify_atom_dict")

    def _subset_atom_dict(self, atom_dict, index_list=None, element_list=None):
        """
        Downselect atom sites to user supplied atom indices or element symbols

        Params
        ------
        atom_dict: Dict
            Dictionary of format {<atom.rdkit.obj>:[atom Symbol, atom Total Valence, atom Degree, # bonded H atoms]}
        index_list: list
            list of integers denoting IUPAC atom indexes of base atoms to be ionized (eg. [0,1])
        element_list: list
            list of atoms symbols denoting base atom types to be ionized (eg. ['O','N'])

        Returns
        -------
        subset_atom_dict defined as {<atom.rdkit.obj>:['C',3,2,1]}
            -atom.GetSymbol() returns atomic symbol e.g. 'C'
            -atom.GetTotalValence() returns total (explicit+implicit) valence
            -atom.GetDegree() returns number of directly-bonded neighbours
             indep. of bond order, dep. of Hs set to explicit
            -atom.GetTotalNumHs returns total (explicit+implicit) Hs on atom

        """

        subset_atom_dict = {}
        if index_list is not None:
            for k, v in atom_dict.items():
                if k in index_list:
                    subset_atom_dict[k] = v

        if element_list is not None:
            for k, v in atom_dict.items():
                if v[0] in element_list:
                    subset_atom_dict[k] = v

        return subset_atom_dict

    def _forcefield_selector(self, forcefield, mw):
        """
        Checks if user specified forcefield is compatible with supplied mol

        Parameters
        ----------
        forcefield: str
            Specify "UFF" for Universal Force Field optimization by RDKit
            Specify "MMFF" or "MMFF94" for Merck Molecular Force Field 94
            Specify "MMFF94s" for the MMFF94 s variant
        mw: RDKit mol object

        Returns
        -------
        RDKit forcefield optimization function that must be implemented

        """

        forcefield = forcefield.lower()
        if forcefield == "uff":
            if Chem.rdForceFieldHelpers.UFFHasAllMoleculeParams(mw) is True:
                return Chem.AllChem.UFFOptimizeMolecule
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

    def _update_geometry_charge(self, geom):
        """ """

        geom.__dict__.update(charge=geom.get_formal_charge())

    def _add_ion(
        self, init_mol, base_atom_idx, ion_atomic_num, sanitize, forcefield, ff_iter
    ):
        """
        Adds specified ion (by atomic num) to specified base atom (by atom's index).

        Parameters
        ----------
        init_mol: RDKit mol object
        base_atom_idx: Int
            Integer denoting the the atom to which an ion will be added
        ion_atomic_num: Int
            Integer denoting the atomic number of the ion to be added
        sanitize: Bool
            Kekulize, check valencies, set aromaticity, conjugation and hybridization
        forcefield: str
            Specify "UFF" for Universal Force Field optimization by RDKit
            Specify "MMFF" or "MMFF94" for Merck Molecular Force Field 94
            Specify "MMFF94s" for the MMFF94 s variant
        ff_iter: Int
            Integer specifying the max iterations for the specified forcefield

        Returns
        _______
        Ionized RDKit mol object

        """

        mw = Chem.RWMol(init_mol)
        atom_idx = mw.AddAtom(Chem.Atom(ion_atomic_num))
        ba = mw.GetAtomWithIdx(base_atom_idx)
        fc = ba.GetFormalCharge()
        mw.AddBond(base_atom_idx, atom_idx, Chem.BondType.SINGLE)
        ba.SetFormalCharge(fc + 1)
        mw.UpdatePropertyCache(strict=False)
        if sanitize is True:
            Chem.SanitizeMol(mw)
        if forcefield:
            tempff = self._forcefield_selector(forcefield, mw)
            if "mmff94s" in forcefield:
                tempff(mw, mmffVariant=forcefield, maxIters=ff_iter)
            else:
                tempff(mw, maxIters=ff_iter)
        geom = isicle.geometry.Geometry(mol=mw, basename=self.geom.get_basename())
        self._update_geometry_charge(geom)
        return geom

    def _apply_add_ion(
        self, init_mol, atom_dict, ion_atomic_num, sanitize, forcefield, ff_iter
    ):
        """
        Iteratively ionize all atoms in a supplied dictionary.

        Parameters
        ----------
        init_mol: RDKit mol object
        atom_dict: Dict
            Dictionary of format {<atom.rdkit.obj>:[atom Symbol, atom Total Valence, atom Degree, # bonded H atoms]}
        ion_atomic_num: Int
            Integer denoting the atomic number of the ion to be added
        sanitize: Bool
            Kekulize, check valencies, set aromaticity, conjugation and hybridization
        forcefield: str
            Specify "UFF" for Universal Force Field optimization by RDKit
            Specify "MMFF" or "MMFF94" for Merck Molecular Force Field 94
            Specify "MMFF94s" for the MMFF94 s variant
        ff_iter: Int
            Integer specifying the max iterations for the specified forcefield

        Returns
        -------
        Dictionary of style {base atom index: ionized RDKit mol object}

        """

        mol_dict = {}
        for key in atom_dict.keys():
            mw = self._add_ion(
                init_mol, key, ion_atomic_num, sanitize, forcefield, ff_iter
            )
            # Format {atom_base_idx:<newmol>}
            mol_dict[key] = mw
        return mol_dict

    def _positive_mode(
        self,
        init_mol,
        ion_atomic_num,
        batch,
        index_list,
        element_list,
        sanitize,
        forcefield,
        ff_iter,
        include_Alkali_ne,
    ):
        """
        Subsets atoms in supplied mol by specified index or elemental constraints, or feasibilty based upon supplied ion atomic number.

        Parameters
        ----------
        init_mol: RDKit mol object
        ion_atomic_num: Int
            Integer denoting the atomic number of the ion to be added
        batch: Bool
            Whether adducts should be formed from all possible heavy atoms
        index_list: list
            list of integers denoting IUPAC atom indexes of base atoms to be ionized (eg. [0,1])
        element_list: list
            list of atoms symbols denoting base atom types to be ionized (eg. ['O','N'])
        sanitize: Bool
            Kekulize, check valencies, set aromaticity, conjugation and hybridization
        forcefield: str
            Specify "UFF" for Universal Force Field optimization by RDKit
            Specify "MMFF" or "MMFF94" for Merck Molecular Force Field 94
            Specify "MMFF94s" for the MMFF94 s variant
        ff_iter: Int
            Integer specifying the max iterations for the specified forcefield
        include_Alkaline_ne: Bool
            If False: removes alkali and alkaline metals from base atoms to be ionized
            If True: Adducts are attempted to be formed at these atoms, in addition to other atoms.

        Returns
        -------
        Dictionary of style {base atom index: ionized RDKit mol object}

        """

        atom_dict = self._modify_atom_dict(
            init_mol,
            ion_atomic_num,
            mode="positive",
            include_Alkali_ne=include_Alkali_ne,
        )
        if index_list is not None:
            subset_atom_dict = self._subset_atom_dict(atom_dict, index_list=index_list)
            mol_dict = self._apply_add_ion(
                init_mol,
                subset_atom_dict,
                ion_atomic_num,
                sanitize,
                forcefield,
                ff_iter,
            )
        elif element_list is not None:
            subset_atom_dict = self._subset_atom_dict(
                atom_dict, element_list=element_list
            )
            mol_dict = self._apply_add_ion(
                init_mol,
                subset_atom_dict,
                ion_atomic_num,
                sanitize,
                forcefield,
                ff_iter,
            )
        elif batch == True:
            mol_dict = self._apply_add_ion(
                init_mol, atom_dict, ion_atomic_num, sanitize, forcefield, ff_iter
            )
        else:
            raise ValueError(
                "Must provide a list of atom indexes, list of elements, or set batch=True."
            )

        return mol_dict

    def _remove_ion(
        self, init_mol, base_atom_idx, ion_atomic_num, sanitize, forcefield, ff_iter
    ):
        """
        Removes specified ion (by atomic num) from specified base atom (by atom's index).

        Parameters
        ----------
        init_mol: RDKit mol object
        base_atom_idx: Int
            Integer denoting the the atom from which an ion will be removed from
        ion_atomic_num: Int
            Integer denoting the atomic number of the ion to be removed
        sanitize: Bool
            Kekulize, check valencies, set aromaticity, conjugation and hybridization
        forcefield: str
            Specify "UFF" for Universal Force Field optimization by RDKit
            Specify "MMFF" or "MMFF94" for Merck Molecular Force Field 94
            Specify "MMFF94s" for the MMFF94 s variant
        ff_iter: Int
            Integer specifying the max iterations for the specified forcefield

        Returns
        -------
        Deionized RDKit mol object

        """

        mw = Chem.RWMol(init_mol)
        ba = mw.GetAtomWithIdx(base_atom_idx)
        fc = ba.GetFormalCharge()
        nbrs_dict = {nbr.GetAtomicNum(): nbr.GetIdx() for nbr in ba.GetNeighbors()}
        mw.RemoveAtom(nbrs_dict[ion_atomic_num])
        ba.SetFormalCharge(fc - 1)
        mw.UpdatePropertyCache(strict=False)
        if sanitize is True:
            Chem.SanitizeMol(mw)
        if forcefield:
            tempff = self._forcefield_selector(forcefield, mw)
            if "mmff94s" in forcefield:
                tempff(mw, mmffVariant=forcefield, maxIters=ff_iter)
            else:
                tempff(mw, maxIters=ff_iter)
        geom = isicle.geometry.Geometry(mol=mw, basename=self.geom.get_basename())
        self._update_geometry_charge(geom)
        return geom

    def _apply_remove_ion(
        self, init_mol, atom_dict, ion_atomic_num, sanitize, forcefield, ff_iter
    ):
        """
        Iteratively deionize all atoms in a supplied dictionary.

        Parameters
        ----------
        init_mol: RDKit mol object
        atom_dict: Dict
            Dictionary of format {<atom.rdkit.obj>:[atom Symbol, atom Total Valence, atom Degree, # bonded H atoms]}
        ion_atomic_num: Int
            Integer denoting the atomic number of the ion to be removed
        sanitize: Bool
            Kekulize, check valencies, set aromaticity, conjugation and hybridization
        forcefield: str
            Specify "UFF" for Universal Force Field optimization by RDKit
            Specify "MMFF" or "MMFF94" for Merck Molecular Force Field 94
            Specify "MMFF94s" for the MMFF94 s variant
        ff_iter: Int
            Integer specifying the max iterations for the specified forcefield

        Returns
        -------
        Dictionary of style {base atom index: deionized RDKit mol object}

        """

        mol_dict = {}
        for key, value in atom_dict.items():
            mw = self._remove_ion(
                init_mol, key, ion_atomic_num, sanitize, forcefield, ff_iter
            )
            mol_dict[key] = mw
        return mol_dict

    def _negative_mode(
        self,
        init_mol,
        ion_atomic_num,
        batch,
        index_list,
        element_list,
        sanitize,
        forcefield,
        ff_iter,
        include_Alkali_ne,
    ):
        """
        Subsets atoms in supplied mol by specified index or elemental constraints, or feasibilty based upon supplied ion atomic number.

        Parameters
        ----------
        init_mol: RDKit mol object
        ion_atomic_num: Int
            Integer denoting the atomic number of the ion to be removed
        batch: Bool
            Whether adducts should be formed from all possible heavy atoms
        index_list: list
            list of integers denoting IUPAC atom indexes of base atoms to be deionized (eg. [0,1])
        element_list: list
            list of atoms symbols denoting base atom types to be ionized (eg. ['O','N'])
        sanitize: Bool
            Kekulize, check valencies, set aromaticity, conjugation and hybridization
        forcefield: str
            Specify "UFF" for Universal Force Field optimization by RDKit
            Specify "MMFF" or "MMFF94" for Merck Molecular Force Field 94
            Specify "MMFF94s" for the MMFF94 s variant
        ff_iter: Int
            Integer specifying the max iterations for the specified forcefield
        include_Alkaline_ne: Bool
            If False: removes alkali and alkaline metals from base atoms to be deionized
            If True: Adducts are attempted to be formed at these atoms, in addition to other atoms.

        Returns
        -------
        Dictionary of style {base atom index: deionized RDKit mol object}

        """

        atom_dict = self._modify_atom_dict(
            init_mol,
            ion_atomic_num,
            mode="negative",
            include_Alkali_ne=include_Alkali_ne,
        )
        if index_list is not None:
            subset_atom_dict = self._subset_atom_dict(atom_dict, index_list=index_list)
            mol_dict = self._apply_remove_ion(
                init_mol,
                subset_atom_dict,
                ion_atomic_num,
                sanitize,
                forcefield,
                ff_iter,
            )
        elif element_list is not None:
            subset_atom_dict = self._subset_atom_dict(
                atom_dict, element_list=element_list
            )
            mol_dict = self._apply_remove_ion(
                init_mol,
                subset_atom_dict,
                ion_atomic_num,
                sanitize,
                forcefield,
                ff_iter,
            )
        elif batch == True:
            mol_dict = self._apply_remove_ion(
                init_mol, atom_dict, ion_atomic_num, sanitize, forcefield, ff_iter
            )
        else:
            raise ValueError(
                "Must provide a list of atom indexes, list of elements, or set batch=True."
            )
        return mol_dict

    def submit(
        self,
        batch=True,
        index_list=None,
        element_list=None,
        sanitize=True,
        forcefield="UFF",
        ff_iter=200,
        include_Alkali_ne=False,
    ):
        """
        Calls positive or negative ionization modes based upon supplied ion list.

        Parameters
        ----------
        batch: Bool
            Whether adducts should be formed from all possible heavy atoms
        index_list: list
            list of integers denoting IUPAC atom indexes of base atoms to be ionized (eg. [0,1])
        element_list: list
            list of atoms symbols denoting base atom types to be ionized (eg. [`O`,`N`])
        sanitize: Bool
            Kekulize, check valencies, set aromaticity, conjugation and hybridization
        forcefield: str
            Specify `UFF` for Universal Force Field optimization by RDKit (default)
            Specify `MMFF` or `MMFF94` for Merck Molecular Force Field 94
            Specify `MMFF94s` for the MMFF94 s variant
        ff_iter: Int
            Integer specifying the max iterations for the specified forcefield (200 default)
        include_Alkaline_ne: Bool
            If False: removes alkali and alkaline metals from base atoms to be ionized (default)
            If True: Adducts are attempted to be formed at these atoms, in addition to other atoms.

        Returns
        -------
        Dictionary of style {base atom index: ionized RDKit mol object}

        """

        pt = Chem.GetPeriodicTable()

        anion_dict = {}
        cation_dict = {}
        complex_dict = {}

        if index_list is not None:
            index_list = safelist(index_list)
        if element_list is not None:
            element_list = safelist(element_list)

        for x in self.adducts["cations"]:
            ion_atomic_num = pt.GetAtomicNumber(re.findall("(.+?)[0-9]?[+]", x)[0])
            mol_dict = self._positive_mode(
                self.geom.mol,
                ion_atomic_num,
                batch,
                index_list,
                element_list,
                sanitize,
                forcefield,
                ff_iter,
                include_Alkali_ne,
            )
            cation_dict[x] = list(mol_dict.values())

        for x in self.adducts["anions"]:
            ion_atomic_num = pt.GetAtomicNumber(re.findall("(.+?)[0-9]?[-]", x)[0])
            mol_dict = self._negative_mode(
                self.geom.mol,
                ion_atomic_num,
                batch,
                index_list,
                element_list,
                sanitize,
                forcefield,
                ff_iter,
                include_Alkali_ne,
            )
            anion_dict[x] = list(mol_dict.values())

        for x in self.adducts["complex"]:
            # This self-determined charge doesn't account for K(2+) or CO2- type adducts
            parsed = re.findall(".*?[+-]", x)
            # TODO strip whitespace
            geom_list = safelist(self.geom)

            for ion in parsed:
                ion_atomic_num = pt.GetAtomicNumber(
                    re.findall("(.+?)[0-9]?[+-]", ion)[0]
                )
                if ("+") in ion:
                    res = list(
                        map(
                            lambda geom: list(
                                self._positive_mode(
                                    geom.mol,
                                    ion_atomic_num,
                                    batch,
                                    index_list,
                                    element_list,
                                    sanitize,
                                    forcefield,
                                    ff_iter,
                                    include_Alkali_ne,
                                ).values()
                            ),
                            geom_list,
                        )
                    )
                    complex_dict[x] = functools.reduce(operator.iconcat, res, [])

                elif ("-") in ion:
                    res = list(
                        map(
                            lambda geom: list(
                                self._negative_mode(
                                    geom.mol,
                                    ion_atomic_num,
                                    batch,
                                    index_list,
                                    element_list,
                                    sanitize,
                                    forcefield,
                                    ff_iter,
                                    include_Alkali_ne,
                                ).values()
                            ),
                            geom_list,
                        )
                    )
                    complex_dict[x] = functools.reduce(operator.iconcat, res, [])
                geom_list = complex_dict.get(x)

        # Combine ion dictionaries and index
        ensemble = []
        for k, v in {**cation_dict, **anion_dict, **complex_dict}.items():
            id = 0
            for geom in v:
                geom.__dict__.update(ion=k, adductID=id)
                ensemble.append(geom)
                id += 1

        self.adducts = build_adduct_ensemble(ensemble)

    def finish(self):
        """
        Parse results.

        Returns
        -------
        :obj:`~isicle.adducts.ExplicitIonizationWrapper`

        """

        return self

    def run(self, geom, ion_path=None, ion_list=None, **kwargs):
        """
        Ionize geometry via with supplied geometry and file containing list of ions.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        ion_path : str
            Filepath to text file containing ions with charge (eg. `H+`) to be considered
            Either ion_path or ion_list must be specified
        ion_list : list
            List of strings of adducts to be considered.
            Must be specifed in syntax `Atom+` or `Atom-`, eg. `H+`, `Na+`, `H-Na+`
            Either ion_path or ion_list must be specified
        **kwargs
            Keyword arguments to configure how ionization is run.
            See :meth:`~isicle.adducts.ExplicitIonizationWrapper.submit`.

        Returns
        -------
        :obj:`~isicle.adducts.ExplicitIonizationWrapper`

        """

        # New instance
        self = ExplicitIonizationWrapper(**geom.__dict__)

        # Sets geom to self.geom
        self.set_geometry(geom)

        # Infers charge of geom.mol
        self.set_charge()

        # Load specified ions by type
        # Validity check if negative ionization can be done
        self.configure(ion_path=ion_path, ion_list=ion_list)

        # Generate adducts
        self.submit(**kwargs)

        # Finish
        self.finish()

        return self

    def get_structures(self):
        """
        Extract all structures from containing object as a conformational ensemble.

        Returns
        -------
        :obj:`~isicle.conformers.ConformationalEnsemble`
            Conformational ensemble.

        """

        return self.adducts.get_structures()

    def get_adducts(self):
        """
        Extract all structures from containing object as an adduct ensemble.

        Returns
        -------
        :obj:`~isicle.adducts.AdductEnsemble`
            Adduct ensemble.

        """

        return self.adducts

    def save_pickle(self, path):
        """
        Pickle this class instance.

        Parameters
        ----------
        path : str
            Path to save pickle object to.

        """

        with open(path, "wb") as f:
            pickle.dump(self, f)

    def save(self, path):
        """
        Save this class instance.

        Parameters
        ----------
        path : str
            Path to save object to.

        """

        self.save_pickle(path)


class CRESTIonizationWrapper(WrapperInterface):
    """ """

    _defaults = ("geom", "adducts")
    _default_value = None

    def __init__(self, **kwargs):
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))
        self.__dict__.update(**kwargs)
        if self.adducts is None:
            self.adducts = {}

    def set_geometry(self, geom):
        """
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        """

        if self.__dict__.get("geom") is not None:
            pass
        elif self.__dict__.get("mol") is not None:
            self.geom = geom
        elif self.__dict__.get("xyz") is not None:
            self.geom = geom
        else:
            raise ValueError("Could not find self.xyz, self.mol, or self.geom")

    def set_charge(self):
        """ """

        if self.geom.__dict__.get("charge") is None:
            if self.geom.load.get("filetype") == ".xyz":
                # self.geom.__dict__.update(charge=charge)
                raise ValueError(
                    "Must first run geom.set_formal_charge for an xyz structure"
                )
            else:
                self.geom.__dict__.update(charge=self.geom.get_formal_charge())

    def _set_ions(self, ion_path=None, ion_list=None):
        cations, anions, complex = _parse_ions(ion_path=ion_path, ion_list=ion_list)
        # cast to list safely
        self.adducts["cations"] = safelist(cations)
        self.adducts["anions"] = safelist(anions)
        self.adducts["complex"] = safelist(complex)

    def _filter_supported_by_xtb(self, unknown_valid_list):
        """
        Filters ions by what is supported by xtb software.
        consult https://github.com/grimme-lab/crest/issues/2 for supported ions by xtb CREST
        xtb supports alkali and alkaline earth metals.

        """

        pt = Chem.GetPeriodicTable()
        valid_list = []
        for ion in unknown_valid_list:
            parsed_ion = re.findall("(.+?)[0-9]?[+-]", ion)
            if len(parsed_ion) == 0:
                raise ValueError(
                    f"Couldn't identify supplied ion {ion} in _filter_supported_by_xtb"
                )
            else:
                parsed_num = [pt.GetAtomicNumber(i) for i in parsed_ion]
                group_check = [_check_atom_group(i) for i in parsed_num]
                if False in group_check:
                    continue
                else:
                    valid_list.append(ion)
        return valid_list

    def _check_valid(self):
        """
        Performs substructure search and checks if supported by CREST documentation

        """

        self.adducts["cations"] = self._filter_supported_by_xtb(self.adducts["cations"])
        self.adducts["anions"] = self._filter_supported_by_xtb(self.adducts["anions"])
        self.adducts["complex"] = self._filter_supported_by_xtb(self.adducts["complex"])

        if self.geom.__dict__.get("mol") is not None:
            self.adducts["anions"] = _filter_by_substructure_match(
                self.geom.mol, self.adducts["anions"]
            )
            self.adducts["complex"] = _filter_by_substructure_match(
                self.geom.mol, self.adducts["complex"]
            )

    def configure(self, ion_path=None, ion_list=None):
        self._set_ions(ion_path=ion_path, ion_list=ion_list)
        self._check_valid()

    def _parse_ion_charge(self, ion):
        """ """

        # TODO update with MSAC code for ion charge lookup table
        charge = re.findall("[0-9]", ion)
        if not charge:
            charge = 1
        else:
            charge = int(charge[0])
        return charge

    def _update_geometry_charge(self, geom, adduct, ion_charge, mode):
        """ """

        if mode == "negative":
            charge = geom.__dict__.get("charge") - ion_charge
        elif mode == "positive":
            charge = geom.__dict__.get("charge") + ion_charge
        adduct.__dict__.update(charge=charge)

    def _positive_mode(
        self,
        geom,
        forcefield,
        ewin,
        cation,
        optlevel,
        dryrun,
        processes,
        solvation,
        ignore_topology,
    ):
        """
        Call isicle.md.md for specified geom and cation
        """

        output = md(
            geom,
            program="xtb",
            task="protonate",
            forcefield=forcefield,
            ewin=ewin,
            ion=cation,
            optlevel=optlevel,
            dryrun=dryrun,
            charge=geom.charge,
            processes=processes,
            solvation=solvation,
            ignore_topology=ignore_topology,
        )
        ion_charge = self._parse_ion_charge(cation)
        for adduct in output.geom:
            self._update_geometry_charge(geom, adduct, ion_charge, "positive")
        return output.geom

    def _negative_mode(
        self,
        geom,
        forcefield,
        ewin,
        anion,
        optlevel,
        dryrun,
        processes,
        solvation,
        ignore_topology,
    ):
        """
        Call isicle.md.md for specified geom and anion
        """

        output = md(
            geom,
            program="xtb",
            task="deprotonate",
            forcefield=forcefield,
            ewin=ewin,
            ion=anion,
            optlevel=optlevel,
            dryrun=dryrun,
            charge=geom.charge,
            processes=processes,
            solvation=solvation,
            ignore_topology=ignore_topology,
        )
        ion_charge = self._parse_ion_charge(anion)
        for adduct in output.geom:
            self._update_geometry_charge(geom, adduct, ion_charge, "negative")
        return output.geom

    def submit(
        self,
        forcefield="gfn2",
        ewin=30,
        optlevel="Normal",
        dryrun=False,
        processes=1,
        solvation=None,
        ignore_topology=False,
    ):
        """
        Call positive_mode and negative_mode to ionize according to parsed ion lists
        isicle.md.md returned by both modes to self.anion, self.cation, self.complex \
        in form of {ion : xyz object [returned by xtb]}

        Input
        -----
        see isicle.md.md for forcefield, ewin, optlevel, dryrun
        molecular_charge is net charge of structure pre-ionization, default neutral, 0

        """

        anion_dict = {}
        cation_dict = {}
        complex_dict = {}

        for x in self.adducts["cations"]:
            cation_dict[x] = self._positive_mode(
                self.geom,
                forcefield,
                ewin,
                x,
                optlevel,
                dryrun,
                processes,
                solvation,
                ignore_topology,
            )

        for x in self.adducts["anions"]:
            anion_dict[x] = self._negative_mode(
                self.geom,
                forcefield,
                ewin,
                x,
                optlevel,
                dryrun,
                processes,
                solvation,
                ignore_topology,
            )

        for x in self.adducts["complex"]:
            # This self-determined charge doesn't account for K(2+) or CO2- type adducts
            parsed = re.findall(".*?[+-]", x)
            # TODO strip whitespace
            geom_list = safelist(self.geom)
            for ion in parsed:
                if ("+") in ion:
                    # complex_dict[x] =
                    res = list(
                        map(
                            lambda geom: self._positive_mode(
                                geom,
                                forcefield,
                                ewin,
                                ion,
                                optlevel,
                                dryrun,
                                processes,
                                solvation,
                                ignore_topology,
                            ),
                            geom_list,
                        )
                    )
                    complex_dict[x] = functools.reduce(operator.iconcat, res, [])

                elif ("-") in ion:
                    res = list(
                        map(
                            lambda geom: self._negative_mode(
                                geom,
                                forcefield,
                                ewin,
                                ion,
                                optlevel,
                                dryrun,
                                processes,
                                solvation,
                                ignore_topology,
                            ),
                            geom_list,
                        )
                    )
                    complex_dict[x] = functools.reduce(operator.iconcat, res, [])
                geom_list = complex_dict.get(x)

        # Combine ion dictionaries and index
        ensemble = []
        for k, v in {**cation_dict, **anion_dict, **complex_dict}.items():
            id = 0
            for geom in v:
                geom.__dict__.update(ion=k, adductID=id)
                ensemble.append(geom)
                id += 1

        self.adducts = build_adduct_ensemble(ensemble)

    def finish(self):
        """
        Parse results.

        Returns
        -------
        :obj:`~isicle.adducts.CRESTIonizationWrapper`

        """

        return self

    def run(self, geom, ion_path=None, ion_list=None, **kwargs):
        """
        Ionize geometry via with supplied geometry and file containing list of ions.
        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        ion_method : str
            Alias for ionaztion method selection (explicit).
        ion_list : list
            List of ion str.
            See :meth`~isicle.adducts.parse_ions`.
        **kwargs
            Keyword arguments to configure how ionization is run.
            See :meth:`~isicle.adducts.CRESTIonizationWrapper.generator`
        Returns
        -------
        :obj:`~isicle.adducts.CRESTIonizationWrapper`

        """

        self = CRESTIonizationWrapper(**geom.__dict__)

        self.set_geometry(geom)

        # Infers charge for non xyz files, infers default neutral for xyz files
        self.set_charge()

        # Load specified ions by type
        # Validity check if negative ionization can be done
        # Validity check if CREST supports the ions specified
        self.configure(ion_path=ion_path, ion_list=ion_list)

        # Generate adducts
        self.submit(**kwargs)

        # Compile various adducts into one dictionary: self.adducts
        self.finish()

        return self

    def get_structures(self):
        """
        Extract all structures from containing object as a conformational ensemble.

        Returns
        -------
        :obj:`~isicle.conformers.ConformationalEnsemble`
            Conformational ensemble.

        """

        return self.adducts.get_structures()

    def get_adducts(self):
        """
        Extract all structures from containing object as an adduct ensemble.

        Returns
        -------
        :obj:`~isicle.adducts.AdductEnsemble`
            Adduct ensemble.

        """

        return self.adducts

    def save_pickle(self, path):
        """
        Pickle this class instance.

        Parameters
        ----------
        path : str
            Path to save pickle object to.

        """

        with open(path, "wb") as f:
            pickle.dump(self, f)

    def save(self, path):
        """
        Save this class instance.

        Parameters
        ----------
        path : str
            Path to save object to.

        """

        self.save_pickle(path)
