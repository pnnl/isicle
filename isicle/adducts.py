import isicle
from isicle.interfaces import WrapperInterface
from isicle.md import md
from isicle.utils import safelist
from isicle.conformers import ConformationalEnsemble
import os
from rdkit import Chem
import re
import rdkit.Chem.rdchem as rdc


def _parse_file_ions(path):
    '''
    Read adduct ions (eg. 'Na+', 'H+', 'Mg2+') from file.

    Parameters
    ----------
    path : str
        Path to ions text file (.txt, etc.), each ion should be a new line.

    Returns
    -------
    parsed_contents
        List of ions from given text file.
    '''
    if path is None:
        raise TypeError('Path to file containing ions must be specified')
    contents = geometry._load_text(path)
    parsed_contents = [line.rstrip() for line in contents]
    return parsed_contents


def _parse_list_ions(ion_list):
    '''
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
    '''
    anions = []
    cations = []
    complex = []
    ion_list = safelist(ion_list)
    for x in ion_list:
        if ('+' in x) and ('-' in x):
            complex.append(x)
        elif '+' in x:
            cations.append(x)
        elif '-' in x:
            anions.append(x)
        else:
            raise ValueError('Unrecognized ion specified: {}'.format(x))
    return (cations, anions, complex)


def _parse_ions(ion_path=None, ion_list=None):
    if ion_path is not None:
        # Load ion file
        ion_list = _parse_file_ions(ion_path)
    if ion_list is not None:
        # Parse ion file
        return _parse_list_ions(ion_list)
    else:
        raise RuntimeError('No ions supplied to parse.')


def ionize(ion_method):
    '''
    Selects a supported ionization method.
    Currently only Explicit ionization has been implemented.
    Parameters
    ----------
    ion_method : str
        Alias for ion method selection (e.g. explicit).
    Returns
    -------
    program
        Wrapped functionality of the selected program. Must implement
        :class:`~isicle.interfaces.AdductInterface`.
    '''
    ion_method_map = {'explicit': ExplicitIonizationWrapper, 'crest': CRESTIonizationWrapper}

    if ion_method.lower() in ion_method_map.keys():
        return ion_method_map[ion_method.lower()]
    else:
        raise ValueError('{} not a supported ionization method.'.format(ion_method))


def _check_atom_group(ion_atomic_num):
    """
    Checks periodic group atom belongs to.

    """
    return ion_atomic_num in [1, 3, 4, 11, 12, 19, 20, 37, 38, 55, 56, 87, 88]


def _filter_by_substructure_match(init_mol, unknown_valid_list):
    valid_list = []
    for ion in unknown_valid_list:
        parsed_ion = re.findall('.+?[0-9]?[+-]', ion)
        subset = [f"[{re.findall('(.+?)[0-9]?[-]', i)[0]}]" for i in parsed_ion if '-' in i]
        subset_mol = [Chem.MolFromSmiles(i, sanitize=False) for i in subset]
        substruc_check = [i for i in subset_mol if init_mol.HasSubstructMatch(i)]
        if False in substruc_check:
            continue
        else:
            valid_list.append(ion)
    return valid_list


def write(IonizationWrapper, path=None, fmt=None):
    '''
    Write mol objects to file.

    Parameters
    ----------
    write_files : boolean
        Indicate whether to write all mol objects to file
    path : str
        Directory to write output files. Only used if `write_files` is
        True
    fmt : str
        Format in which to save the RDKit mol object
    '''
    if (path is not None) and (fmt is not None):
        if path is None:
            raise ValueError('Must supply `path`.')
        if fmt is None:
            raise ValueError('Must supply `fmt`.')
        for key, value in IonizatinWrapper.adducts.items():
            for key2, value2 in value.items():
                isicle.geometry.Geometry.save(value2, os.path.join(
                    path, '{}{}.{}'.format(key, key2, fmt)), fmt)
        return
    elif path is not None:
        raise(
            'path passed to isicle.adducts.ExplicitIonizationWrapper.finish; fmt flag must also be passed; data not saved.')
    elif fmt is not None:
        raise(
            'fmt is supplied to isicle.adducts.ExplicitIonizationWrapper.finish; path flag must also be passed; data not saved. ')
    else:
        return


def proton_affinity(MH, M, temp=298.15):
    '''
    Calculate proton affinity relative to passed M+H, M values
    𝑃𝐴(𝐵)=[𝐸𝑒𝑙𝑒(𝐵)−𝐸𝑒𝑙𝑒(𝐵𝐻+)]+[𝑍𝑃𝐸(𝐵)−𝑍𝑃𝐸(𝐵𝐻+)]+(5/2)𝑅𝑇

    Parameters
    ----------
    MH : dict
        Dict with keys: energy (kcal/mol), zpe (kcal/mol)
    M : dict
        Dict with keys(): energy (kcal/mol), zpe (kcal/mol)
    temp : opt, default 298.15 K
    '''
    R = 0.00198720425864083
    PA = M['energy'] - MH['energy'] + M['zpe'] - MH['zpe'] + (5/2)*R*temp
    return PA


def gas_basicity(MH, M, temp=298.15, SH=108.8):
    '''
    Calculate gas basicity relative to passed M+H, M
    𝐺𝐵(𝐵)=[𝐸𝑒𝑙𝑒(𝐵)−𝐸𝑒𝑙𝑒(𝐵𝐻+)]+[𝑍𝑃𝐸(𝐵)−𝑍𝑃𝐸(𝐵𝐻+)]+(5/2)𝑅𝑇−𝑇[𝑆(𝐵)+𝑆(𝐻+)−𝑆(𝐵𝐻+)]

    Parameters
    ----------
    MH : dict
        Dict with keys: zpe, enthalpy, entropy
    M : dict
        Dict with keys(): zpe, enthalpy, entropy
    temp : opt , default 298.15 K
    SH : opt, default 108.8 J.mol/K
    '''
    # need Energy, ZPE, Temperature, entropy
    # S(H+) = 108.8 J mol/K, 298 K)
    R = 0.00198720425864083
    GB = M['energy'] - MH['energy'] + M['zpe'] - MH['zpe'] + \
        (5/2)*R*temp - temp*(M['entropy'] + SH - MH['entropy'])
    return GB


def build_adduct_ensembl(geometries):
    '''
    Create an adduct ensemble from a collection of geometries.
    Parses adduct dictionaries and returns a list of all geometry objects.

    Parameters
    ----------
    geometries : dictionary object from generator output
                 {'IonCharge' : [single_geom...]}

    Returns
    -------
    list : geometries nested inside adduct dictionaries

    '''
    ensembl = []
    for k, v in geometries.items():
        for geom in v:
            ensembl.append(geom.mol)

    return ensembl


class ExplicitIonizationWrapper(WrapperInterface):
    def __init__(self):
        '''
        Initialize :obj:`~isicle.adducts.ExplicitIonizationWrapper` instance.
        '''
        pass

    def set_geometry(self, geom):
        '''
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.
        Parameters
        ----------
        geometry : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        Note
        ----
        geometry is not modified in isicle.adducts; see
        '''
        # Assign geometry
        self._geom = geom

    def _set_ions(self, ion_path=None, ion_list=None):
        cations, anions, complex = _parse_ions(ion_path=ion_path, ion_list=ion_list)
        # Cast to list safely
        self._cations = safelist(cations)
        self._anions = safelist(anions)
        self._complex = safelist(complex)

    def _check_valid(self):
        '''
        Perform substructure search
        '''
        self._anions = _filter_by_substructure_match(self._geom.mol, self._anions)
        self._complex = _filter_by_substructure_match(self._geom.mol, self._complex)

    def configure(self, ion_path=None, ion_list=None):
        self._set_ions(ion_path=ion_path, ion_list=ion_list)
        self._check_valid()

    def _modify_atom_dict(self, init_mol, ion_atomic_num, mode=None, include_Alkali_ne=False):
        '''
        Downselect atom sites that can accept ions

        Methods
        -------
        all_atom_dict defined as {<atom.rdkit.obj>:['C',3,2,1]}
            -atom.GetSymbol() returns atomic symbol e.g. 'C'
            -atom.GetTotalValence() returns total (explicit+implicit) valence
            -atom.GetDegree() returns number of directly-bonded neighbours
             indep. of bond order, dep. of Hs set to explicit
            -atom.GetTotalNumHs returns total (explicit+implicit) Hs on atom
        '''
        pt = Chem.GetPeriodicTable()
        all_atom_dict = {atom.GetIdx(): [atom.GetSymbol(), atom.GetTotalValence() -
                                         atom.GetDegree(), atom.GetTotalNumHs(includeNeighbors=True),
                                         atom] for atom in init_mol.GetAtoms()}

        if mode == 'positive':
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

        elif mode == 'negative':
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
            raise NotImplementedError

    def _add_ion(self, init_mol, base_atom_idx, ion_atomic_num, sanitize, optimize):
        '''
        Adds specified ion (by atomic num) to specified base atom (by atom's index).
        Returns ionized RDKit mol object.
        '''
        mw = Chem.RWMol(init_mol)
        atom_idx = mw.AddAtom(Chem.Atom(ion_atomic_num))
        ba = mw.GetAtomWithIdx(base_atom_idx)
        fc = ba.GetFormalCharge()
        mw.AddBond(base_atom_idx, atom_idx, Chem.BondType.SINGLE)
        ba.SetFormalCharge(fc+1)
        mw.UpdatePropertyCache(strict=False)
        if sanitize is True:
            Chem.SanitizeMol(mw)
        if optimize is True:
            # TODO revisit EmbedMolecule error for large molecules
            Chem.AllChem.EmbedMolecule(mw)
            Chem.AllChem.UFFOptimizeMolecule(mw)
        geom = isicle.geometry.Geometry(mol=mw, history=self._geom.get_history())
        geom._update_history('ionize')
        return geom

    def _apply_add_ion(self, init_mol, atom_dict, ion_atomic_num, sanitize, optimize):
        '''
        Iterates through base atoms of RDKit mol object.
        Calls method to ionize and return new RDKit mol object.
        Returns dictionary of {base atom: ionized RDKit mol object}.
        '''

        mol_dict = {}
        for key in atom_dict.keys():
            mw = self._add_ion(init_mol, key, ion_atomic_num, sanitize, optimize)
            # Format {atom_base_idx:<newmol>}
            mol_dict[key] = mw
        return mol_dict

    def _positive_mode(self, init_mol, ion_atomic_num, batch, single_atom_idx, sanitize, optimize, include_Alkali_ne):
        '''
        Adds specified cation to RDKit mol object.
        batch=True: adduct generated at every potential base atom
        batch=False, single_atom_idx!=None: adduct generated at specified atom index
        '''
        atom_dict = self._modify_atom_dict(
            init_mol, ion_atomic_num, mode='positive', include_Alkali_ne=include_Alkali_ne)

        if ((batch == False) and (single_atom_idx is not None)):
            mol_dict = self._apply_add_ion(
                init_mol, {single_atom_idx: atom_dict[single_atom_idx]}, ion_atomic_num, sanitize, optimize)
        elif batch == True:
            mol_dict = self._apply_add_ion(init_mol, atom_dict, ion_atomic_num, sanitize, optimize)
        else:
            raise ValueError('Must provide an index value with batch==False.')
        return mol_dict

    def _remove_ion(self, init_mol, base_atom_idx, ion_atomic_num, sanitize, optimize):
        '''
        Removes specified ion (by atomic num) from specified base atom (by atom's index).
        Returns deionized RDKit mol object.
        '''

        mw = Chem.RWMol(init_mol)
        ba = mw.GetAtomWithIdx(base_atom_idx)
        fc = ba.GetFormalCharge()
        nbrs_dict = {nbr.GetAtomicNum(): nbr.GetIdx() for nbr in ba.GetNeighbors()}
        mw.RemoveAtom(nbrs_dict[ion_atomic_num])
        ba.SetFormalCharge(fc-1)
        mw.UpdatePropertyCache(strict=False)
        if sanitize is True:
            Chem.SanitizeMol(mw)
        if optimize is True:
            # TODO revisit EmbedMolecule issues for large molecules
            Chem.AllChem.EmbedMolecule(mw)
            Chem.AllChem.UFFOptimizeMolecule(mw)
        geom = isicle.geometry.Geometry(mol=mw, history=self._geom.get_history())
        geom._update_history('ionize')
        return geom

    def _apply_remove_ion(self, init_mol, atom_dict, ion_atomic_num, sanitize, optimize):
        '''
        Iterates through base atoms of RDKit mol object.
        Calls method to deionize and return new RDKit mol object.
        Returns dictionary of {base atom: deionized RDKit mol object}.
        '''

        mol_dict = {}
        for key, value in atom_dict.items():
            mw = self._remove_ion(
                init_mol, key, ion_atomic_num, sanitize, optimize)
            mol_dict[key] = mw
        return mol_dict

    def _negative_mode(self, init_mol, ion_atomic_num, batch, single_atom_idx, sanitize, optimize, include_Alkali_ne):
        '''
        Removes specified anion from RDKit mol object at every potential base atom.
        '''
        atom_dict = self._modify_atom_dict(
            init_mol, ion_atomic_num, mode='negative', include_Alkali_ne=include_Alkali_ne)
        if ((batch == False) and (single_atom_idx is not None)):
            mol_dict = self._apply_remove_ion(
                init_mol, {single_atom_idx: atom_dict[single_atom_idx]}, ion_atomic_num, sanitize, optimize)
        elif batch == True:
            mol_dict = self._apply_remove_ion(
                init_mol, atom_dict, ion_atomic_num, sanitize, optimize)
        else:
            raise ValueError('Must provide an index value with batch==False.')
        return mol_dict

    def submit(self, batch=True, single_atom_idx=None, sanitize=True, optimize=True, include_Alkali_ne=False):
        '''
        Creates specified adducts at every atom site.
        '''
        pt = Chem.GetPeriodicTable()
        anion_dict = {}
        cation_dict = {}
        complex_dict = {}

        for ion in self._cations:
            ion_atomic_num = pt.GetAtomicNumber(re.findall('(.+?)[0-9]?[+]', ion)[0])
            mol_dict = self._positive_mode(self._geom.mol, ion_atomic_num,
                                           batch, single_atom_idx, sanitize, optimize, include_Alkali_ne)
            cation_dict[ion] = list(mol_dict.values())
            # cation_dict defined as {ion+:{base_atom_index: mol}}
        self._cations = cation_dict
        for ion in self._anions:
            ion_atomic_num = pt.GetAtomicNumber(re.findall('(.+?)[0-9]?[-]', ion)[0])
            mol_dict = self._negative_mode(self._geom.mol, ion_atomic_num,
                                           batch, single_atom_idx, sanitize, optimize, include_Alkali_ne)
            anion_dict[ion] = list(mol_dict.values())
            # anion_dict defined as {ion-:{base_atom_index: mol}}
        self._anions = anion_dict

        for ion in self._complex:
            split_ion = re.findall('.+?[0-9]?[+-]', ion)
            temp_mols = [self._geom]
            out_mols = []
            for sion in split_ion:
                ion_atomic_num = pt.GetAtomicNumber(re.findall('(.+?)[0-9]?[+-]', sion)[0])
                for temp_mol in temp_mols:
                    if '+' in sion:
                        mol_dict = self._positive_mode(temp_mol.mol, ion_atomic_num,
                                                       batch, single_atom_idx, sanitize, optimize, include_Alkali_ne)
                    elif '-' in sion:
                        mol_dict = self._negative_mode(temp_mol.mol, ion_atomic_num,
                                                       batch, single_atom_idx, sanitize, optimize, include_Alkali_ne)
                    else:
                        raise ValueError
                    if len(temp_mols) > 1:
                        out_mols.extend(mol_dict.values())
                    elif len(temp_mols) == 1:
                        temp_mols = mol_dict.values()
            complex_dict[ion] = out_mols
        self._complex = complex_dict

    def finish(self):
        '''
        Combine generated structures and optionally save

        Returns
        -------
        :obj:`~isicle.parse.XTBResult`
            Parsed result data.
        '''
        # ion dict format {ion<charge>:mol}
        self.adducts = {**self._cations, **self._anions, **self._complex}

    def run(self, geom, ion_path=None, ion_list=None, **kwargs):
        '''
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
            See :meth:`~isicle.adducts.ExplicitIonizationWrapper.generator`.
        Returns
        -------
        :obj:`~isicle.adducts.ExplicitIonizationWrapper`

        '''
        # New instance
        self = ExplicitIonizationWrapper()

        self.set_geometry(geom)

        # Load specified ions by type
        # Validity check if negative ionization can be done
        # Validity check if CREST supports the ions specified
        self.configure(ion_path=ion_path, ion_list=ion_list)

        # Generate adducts
        self.submit(**kwargs)

        # Combine varions ion mols into one output dictionary
        self.finish()

        return self


class CRESTIonizationWrapper(WrapperInterface):
    def __init__(self):
        pass

    def set_geometry(self, geom):
        '''
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.
        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        '''
        # Assign geometry
        self._geom = geom

    def _set_ions(self, ion_path=None, ion_list=None):
        cations, anions, complex = _parse_ions(ion_path=ion_path, ion_list=ion_list)
        # Cast to list safely
        self._cations = safelist(cations)
        self._anions = safelist(anions)
        self._complex = safelist(complex)

    def _filter_supported_by_xtb(self, unknown_valid_list):
        """
        Filters ions by what is supported by xtb software.
        consult https://github.com/grimme-lab/crest/issues/2 for supported ions by xtb CREST
        xtb supports alkali and alkaline earth metals.
        """
        pt = Chem.GetPeriodicTable()
        valid_list = []
        for ion in unknown_valid_list:
            parsed_ion = re.findall('(.+?)[0-9]?[+-]', ion)
            if len(parsed_ion) == 0:
                raise ValueError(
                    f"Couldn't identify supplied ion {ion} in _filter_supported_by_xtb")
            else:
                parsed_num = [pt.GetAtomicNumber(i) for i in parsed_ion]
                group_check = [_check_atom_group(i) for i in parsed_num]
                if False in group_check:
                    continue
                else:
                    valid_list.append(ion)
        return valid_list

    def _check_valid(self):
        '''
        Performs substructure search and checks if supported by CREST documentation
        '''
        self._cations = self._filter_supported_by_xtb(self._cations)
        self._anions = self._filter_supported_by_xtb(self._anions)
        self._complex = self._filter_supported_by_xtb(self._complex)

        self._anions = _filter_by_substructure_match(self._geom.mol, self._anions)
        self._complex = _filter_by_substructure_match(self._geom.mol, self._complex)

    def configure(self, ion_path=None, ion_list=None):
        self._set_ions(ion_path=ion_path, ion_list=ion_list)
        self._check_valid()

    def _positive_mode(self, forcefield, ewin, cation, charge, optlevel, dryrun, processes, solvation):
        '''
        Call isicle.md.md for specified geom and cation
        '''
        try:
            output = md(self._geom, program='xtb', task='protonate', forcefield=forcefield,
                        ewin=ewin, ion=cation, optlevel=optlevel, dryrun=dryrun, charge=charge, processes=processes, solvation=solvation).result['geom']
        except:
            output = None
        return output

    def _negative_mode(self, forcefield, ewin, anion, charge, optlevel, dryrun, processes, solvation):
        '''
        Call isicle.md.md for specified geom and anion
        '''
        try:
            output = md(self._geom, program='xtb', task='deprotonate', forcefield=forcefield,
                        ewin=ewin, ion=anion, optlevel=optlevel, dryrun=dryrun, charge=charge, processes=processes, solvation=solvation).result['geom']
        except:
            output = None
        return output

    def submit(self, forcefield='gff', ewin=30, charge=0, optlevel='Normal', dryrun=False, processes=1, solvation=None):
        '''
        Call positive_mode and negative_mode to ionize according to parsed ion lists
        isicle.md.md returned by both modes to self.anion, self.cation, self.complex \
        in form of {ion : xyz object [returned by xtb]}

        Input
        -----
        see isicle.md.md for forcefield, ewin, optlevel, dryrun
        molecular_charge is net charge of structure pre-ionization, default neutral, 0
        '''
        anion_dict = {}
        cation_dict = {}
        complex_dict = {}

        for x in self._cations:
            # TODO update how MD is parsed, ensure is geometry object in memory
            cation_dict[x] = self._positive_mode(
                forcefield, ewin, x, charge, optlevel, dryrun, processes, solvation)
        self._cations = cation_dict

        for x in self._anions:
            anion_dict[x] = self._negative_mode(
                forcefield, ewin, x, charge, optlevel, dryrun, processes, solvation)
        self._anions = anion_dict

        for x in self._complex:
            mol_chrg = charge
            # This self-determined charge doesn't account for K(2+) or CO2- type adducts
            parsed = re.findall('.*?[+-]', x)
            # TODO strip whitespace
            for ion in parsed:
                chrg = re.findall('[0-9]', ion)
                if not chrg:
                    chrg = 1
                else:
                    chrg = int(chrg[0])
                if ('+') in ion:
                    complex_dict[x] = self._positive_mode(
                        forcefield, ewin, ion, mol_chrg, optlevel, dryrun, processes, solvation)
                    mol_chrg += chrg
                elif ('-') in ion:
                    complex_dict[x] = self._negative_mode(
                        forcefield, ewin, ion, mol_chrg, optlevel, dryrun, processes, solvation)
                    mol_chrg -= chrg

        self._complex = complex_dict

    def finish(self):
        '''
        Combine generated structures and optionally save

        Returns
        -------
        :obj:`~isicle.parse.XTBResult`
            Parsed result data.
        '''
        # ion dict format {ion<charge>:mol}
        self.adducts = {**self._cations, **self._anions, **self._complex}

    def run(self, geom, ion_path=None, ion_list=None, **kwargs):
        '''
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

        '''
        self = CRESTIonizationWrapper()

        self.set_geometry(geom)

        # Load specified ions by type
        # Validity check if negative ionization can be done
        # Validity check if CREST supports the ions specified
        self.configure(ion_path=ion_path, ion_list=ion_list)

        # Generate adducts
        self.submit(**kwargs)

        # Compile various adducts into one dictionary: self.adducts
        self.finish()

        return self
