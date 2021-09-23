from isicle.interfaces import IonizeWrapperInterface
from isicle import geometry
from isicle.md import md
from isicle.utils import safelist
from isicle.conformers import ConformationalEnsemble
import os
from rdkit import Chem
import re


def load_ions(path):
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


def parse_ions(ion_list):
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
        if ('+' and '-') in x:
            complex.append(x)
        elif ('+') in x:
            cations.append(x)
        elif ('-') in x:
            anions.append(x)
        else:
            raise ValueError('Unrecognized ion specified: {}'.format(x))
            # TODO insert more descriptive error
    return (cations, anions, complex)


def _ionize_method_selector(ion_method):
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
        return ion_method_map[ion_method.lower()]()
    else:
        raise ValueError('{} not a supported ionization method.'.format(ion_method))


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


def build_adduct_ensemble(geometries):
    '''
    Create an adduct ensemble from a collection of geometries.

    Parameters
    ----------
    geometries : dict of :obj:`~isicle.geometry.Geometry` or related subclass
        Collection of geometry instances.

    Returns
    -------
    :obj:`~isicle.adducts.AdductEnsemble`
        Adduct ensemble.

    '''

    return ConformationalEnsemble(geometries)


def ionize(geom, ion_path=None, ion_list=None, ion_method='explicit', **kwargs):
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
         or :meth:`~isicle.adducts.ExplicitIonizationWrapper.generator`.
    Returns
    -------
    :obj:`~isicle.adducts.ExplicitIonizationWrapper`
        Dictionary of
    '''
    if ion_path is not None:
        # Load ion file
        ion_list = load_ions(ion_path)
    if ion_list is not None:
        # Parse ion file
        cations, anions, complex = parse_ions(ion_list)
    else:
        raise RuntimeError('No ions to parse.')

    # Select ionization method
    iw = _ionize_method_selector(ion_method)

    # TODO make sure removeHs=False
    # Load in geometry of geom object
    iw.set_geometry(geom)

    # Load specified ions by type
    iw.set_ions(cations, anions, complex)

    # TODO substructure check to validate specified ion support
    # iw.check_valid()

    # Generate adducts
    iw.generator(**kwargs)

    # Combine/ optionally write files
    res = iw.finish()
    """
    geom = geom._update_structure(False, mol=res.)

    # Erase old properties and add new event and DFT properties
    geom.global_properties = {}
    geom._update_history('md')
    geom = geom.add_global_properties(res.to_dict())
    """
    return res  # geom, res


class ExplicitIonizationWrapper(IonizeWrapperInterface):
    def __init__(self, geom=None, anions=None, cations=None, complex=None):
        '''
        Initialize :obj:`~isicle.adducts.ExplicitIonizationWrapper` instance.
        Creates temporary directory for intermediate files, establishes aliases
        for preconfigured tasks.
        '''
        self.geom = geom
        self.anions = anions
        self.cations = cations
        self.complex = complex

    def set_geometry(self, geom):
        '''
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.
        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        '''
        # Assign geometry
        self.geom = geom

    def set_ions(self, cations, anions, complex):
        # Cast to list safely
        self.cations = safelist(cations)
        self.anions = safelist(anions)
        self.complex = safelist(complex)

    def check_valid(self):
        '''
        Perform substructure search
        Check ion list against CREST nomenclature for ions, supported ion type
        '''
        # consider calc_new_mz from
        # https://github.com/pnnl/mame/blob/24f9fcb19639d5a1c2ca0f788c55d6a2efe18ca6/mame/mame_utils.py
        # consider add_formala, remove_formula from
        # https://github.com/pnnl/mame/blob/24f9fcb19639d5a1c2ca0f788c55d6a2efe18ca6/mame/formula_module.py
        # mol.HasSubstructMatch((ion.strip('+')).strip('-'))
        # complex structure checking, spectre?
        # consider how substructure would break off, not all atoms may be from same region

        # TODO (1) check if ions are in supported list of ions
        # TODO (2) modify self.anions, self.cations, self.complex to conform to supported ions
        return

    def modify_atom_dict(self, init_mol, ion_atomic_num, mode=None):
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
                if atom_atomic_num == ion_atomic_num:
                    # TODO modify logic here
                    # TODO add check for more ion types
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

        # TODO where is this referenced
        elif mode == None:
            raise NotImplementedError

    def _add_ion(self, init_mol, base_atom_idx, ion_atomic_num, sanitize, optimize):
        '''
        Adds specified ion (by atomic num) to specified base atom (by atom's index).
        Returns ionized RDKit mol object.
        '''
        import rdkit.Chem.rdchem as rdc
        from rdkit.Chem import AllChem
        # TODO only works for single atom ion

        mw = rdc.RWMol(init_mol)
        atom_idx = rdc.RWMol.AddAtom(mw, rdc.Atom(ion_atomic_num))
        # atom idx starts at 0; GetNumAtoms returns last idx val+1
        #atom_idx = init_mol.GetNumAtoms()
        ba = mw.GetAtomWithIdx(base_atom_idx)
        fc = ba.GetFormalCharge()
        mw.AddBond(base_atom_idx, atom_idx, Chem.BondType.SINGLE)
        ba.SetFormalCharge(fc+1)
        mw.UpdatePropertyCache(strict=False)
        if sanitize is True:
            Chem.SanitizeMol(mw)
        if optimize is True:
            AllChem.EmbedMolecule(mw)
            AllChem.UFFOptimizeMolecule(mw)
        return mw

    def add_ion(self, init_mol, atom_dict, ion_atomic_num, sanitize, optimize):
        '''
        Iterates through base atoms of RDKit mol object.
        Calls method to ionize and return new RDKit mol object.
        Returns dictionary of {base atom: ionized RDKit mol object}.
        '''

        mol_dict = {}
        for key in atom_dict.keys():
            mw = self._add_ion(init_mol, key, ion_atomic_num, sanitize=sanitize, optimize=optimize)
            # Format {atom_base_idx:<newmol>}
            mol_dict[key] = mw
        return mol_dict

    def positive_mode(self, init_mol, ion_atomic_num, batch=True, single_atom_idx=None, sanitize=False, optimize=False):
        '''
        Adds specified cation to RDKit mol object.
        batch=True: adduct generated at every potential base atom
        batch=False, single_atom_idx!=None: adduct generated at specified atom index
        '''
        atom_dict = self.modify_atom_dict(init_mol, ion_atomic_num, mode='positive')
        if ((batch == False) and (single_atom_idx is not None)):
            mol_dict = self.add_ion(
                init_mol, {single_atom_idx: atom_dict[single_atom_idx]}, ion_atomic_num, sanitize, optimize)
        elif batch == True:
            mol_dict = self.add_ion(init_mol, atom_dict, ion_atomic_num, sanitize, optimize)
        else:
            raise ValueError('Must provide an index value with batch==False.')
        return mol_dict

    def _remove_ion(init_mol, base_atom_idx, ion_atomic_num, sanitize, optimize):
        '''
        Removes specified ion (by atomic num) from specified base atom (by atom's index).
        Returns deionized RDKit mol object.
        '''
        import rdkit.Chem.rdchem as rdc
        mw = rdc.RWMol(init_mol)
        ba = mw.GetAtomWithIdx(base_atom_idx)
        fc = ba.GetFormalCharge()
        nbrs_dict = {nbr.GetAtomicNum(): nbr.GetIdx() for nbr in base_atom.GetNeighbors()}
        mw.RemoveAtom(nbrs_dict[ion_atomic_num])
        ba.SetFormalCharge(fc-1)
        mw.UpdatePropertyCache(strict=False)
        if sanitize is True:
            Chem.SanitizeMol(mw)
        if optimize is True:
            AllChem.EmbedMolecule(mw)
            AllChem.UFFOptimizeMolecule(mw)
        return mw

    def remove_ion(init_mol, atom_dict, ion_atomic_num, sanitize, optimize):
        '''
        Iterates through base atoms of RDKit mol object.
        Calls method to deionize and return new RDKit mol object.
        Returns dictionary of {base atom: deionized RDKit mol object}.
        '''

        mol_dict = {}
        for key, value in atom_dict.items():
            mw = self._remove_ion(init_mol, key, value[3], ion_atomic_num, sanitize, optimize)
            mol_dict[value[3]] = mw
        return mol_dict

    def negative_mode(self, init_mol, ion_atomic_num, batch=True, single_atom_idx=None, sanitize=False, optimize=False):
        '''
        Removes specified anion from RDKit mol object at every potential base atom.
        '''
        atom_dict = self.modify_atom_dict(init_mol, mode='negative')
        if ((batch == False) and (single_atom_idx is not None)):
            mol_dict = self.remove_ion(
                init_mol, {single_atom_idx: atom_dict[single_atom_idx]}, ion_atomic_num, sanitize, optimize)
        elif batch == True:
            mol_dict = self.remove_ion(init_mol, atom_dict, ion_atomic_num, sanitize, optimize)
        else:
            raise ValueError('Must provide an index value with batch==False.')
        return mol_dict

    def generator(self, sanitize=False, optimize=False):
        '''
        Creates specified adducts at every atom site.
        '''
        pt = Chem.GetPeriodicTable()
        anion_dict = {}
        cation_dict = {}
        complex_dict = {}

        for ion in self.cations:
            sion = ion.strip('+')
            ion_atomic_num = pt.GetAtomicNumber(sion)
            mol_dict = self.positive_mode(self.geom.mol, ion_atomic_num,
                                          sanitize=sanitize, optimize=optimize)
            cation_dict[ion] = mol_dict
            # cation_dict defined as {ion+:{base_atom_index: mol}}
        self.cations = cation_dict
        for ion in self.anions:
            sion = ion.strip('-')
            ion_atomic_num = pt.GetAtomicNumber(sion)
            mol_dict = self.negative_mode(self.geom.mol, ion_atomic_num,
                                          sanitize=sanitize, optimize=optimize)
            anion_dict[ion] = mol_dict
            # anion_dict defined as {ion-:{base_atom_index: self.geom}}
            # TODO make change from mol object returned to self.geom
        self.anions = anion_dict

        for ion in self.complex:
            # TODO iterative method for adduct generation
            # Call both positve_mode and negative_mode
            raise NotImplementedError
        self.complex = complex_dict

    def finish(self, write_files=False, path=None, fmt=None):
        '''
        Combine generated structures and optionally save

        Parameters
        ----------
        write_files : boolean
            Indicate whether to write all mol objects to file
        path : str
            Directory to write output files. Only used if `write_files` is
            True
        fmt : str
            Format in which to save the RDKit mol object
        Returns
        -------
        :obj:`~isicle.parse.NWChemResult`
            Parsed result data
        '''
        ion_dict = {**self.cations, **self.anions, **self.complex}
        # ion dict format {ion<charge>:{base_atom_index: mol}}

        if write_files is True:
            if path is None:
                raise ValueError('Must supply `path`.')
            if fmt is None:
                raise ValueError('Must supply `fmt`.')
            for key, value in ion_dict.items():
                for key2, value2 in value.items():
                    geometry.Geometry.save(value2, os.path.join(
                        path, '{}{}.{}'.format(key, key2, fmt)), fmt)
        return ion_dict


class CRESTIonizationWrapper(IonizeWrapperInterface):
    def __init__(self, geom=None, path=None, cations=None, anions=None, complex=None):
        self.geom = geom
        self.path = path
        self.cations = cations
        self.anions = anions
        self.complex = complex

    def set_geometry(self, geom):
        '''
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.
        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        '''
        # Assign geometry
        self.geom = geom

    def set_ions(self, cations, anions, complex):
        # Cast to list safely
        self.cations = safelist(cations)
        self.anions = safelist(anions)
        self.complex = safelist(complex)

    def check_valid(self, ion):
        '''
        Performs substructure search and checks if supported by CREST documentation
        '''
        # mol.HasSubstructMatch((ion.strip('+')).strip('-'))

        # TODO institute ion list check against list of supported CREST ions
        # downfilter by supported ions
        # consult https://github.com/grimme-lab/crest/issues/2 for supported ions by xtb CREST

    def positive_mode(self, forcefield, ewin, cation, molecular_charge, optlevel, dryrun):
        '''
        Call isicle.md.md for specified geom and cation
        '''
        return md(self.geom, program='xtb', task='protonate', forcefield=forcefield,
                  ewin=ewin, ion=cation, optlevel=optlevel, dryrun=dryrun, charge=molecular_charge)

    def negative_mode(self, forcefield, ewin, anion, molecular_charge, optlevel, dryrun):
        '''
        Call isicle.md.md for specified geom and anion
        '''
        return md(self.geom, program='xtb', task='deprotonate', forcefield=forcefield,
                  ewin=ewin, ion=anion, optlevel=optlevel, dryrun=dryrun, charge=molecular_charge)

    def generator(self, molecular_charge=0, forcefield='gff', ewin=100, optlevel='Normal', dryrun=False):
        '''
        Call positive_mode and negative_mode to ionize according to parsed ion lists
        isicle.md.md returned by both modes to self.anion, self.cation, self.complex \
        in form of {ion : xyz object [returned by xtb]}

        Input
        -----
        see isicle.md.md for forcefield, ewin, optlevel, dryrun
        molecular_charge is net charge of structure pre-ionization, default neutral, 0
        '''
        anions = self.anions
        cations = self.cations
        complex = self.complex
        anion_dict = {}
        cation_dict = {}
        complex_dict = {}

        for x in cations:
            # TODO update how MD is parsed, ensure is geometry object in memory
            cation_dict[x] = self.positive_mode(
                forcefield, ewin, x, molecular_charge, optlevel, dryrun)
        self.cations = cation_dict
        for x in anions:
            anion_dict[x] = self.negative_mode(
                forcefield, ewin, x, molecular_charge, optlevel, dryrun)
        self.anions = anion_dict

        for x in complex:
            mol_chrg = molecular_charge
            parsed = re.findall('.*?[+-]', x)
            # TODO strip whitespace
            for ion in parsed:
                chrg = re.findall('[0-9]', ion)
                if not chrg:
                    chrg = 1
                else:
                    chrg = int(chrg[0])
                if ('+') in ion:
                    complex_dict[x] = self.positive_mode(
                        forcefield, ewin, ion, mol_chrg, optlevel, dryrun)
                    mol_chrg += chrg
                elif ('-') in ion:
                    complex_dict[x] = self.negative_mode(
                        forcefield, ewin, ion, mol_chrg, optlevel, dryrun)
                    mol_chrg -= chrg

        self.complex = complex_dict

    def finish(self):
        '''
        Combine generated structures and optionally save

        Returns
        -------
        :obj:`~isicle.parse.NWChemResult`
            Parsed result data.
        '''
        # TODO parse to identify protonated atom location
        ion_dict = {**self.cations, **self.anions, **self.complex}
        # ion dict format {ion<charge>:mol}
        return ion_dict
