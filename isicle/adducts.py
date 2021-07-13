from isicle.interfaces import IonizeInterface
from isicle import geometry
from isicle.md import md
import os
from rdkit import Chem


def load_ions(path: str):
    '''
    Read adduct ions (eg. Na+, H+) from file.

    Parameters
    ----------
    path : str
        Path to ions text file (.txt, etc.), each ion should be a new line.

    Returns
    -------
    parsed_contents
        List of ions from given text file.
    '''
    contents = geometry._load_text(path)
    parsed_contents = [line.rstrip() for line in contents]
    return parsed_contents


def parse_ions(ion_list):
    '''
    Parse and categorize list of ions by ion type.

    Parameters
    ----------
    anions : list
        List of anions, default None
    cations : list
        List of cations, default None
    complex : list
        List of complex ions, default None

    Returns
    -------

    '''
    anions = []
    cations = []
    complex = []

    # TODO make parser more robust for input list variation
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
        raise ValueError(('{} not a supported ionization method.').format(ion_method))


def ionize(geom, ion_path=None, ion_method='explicit', **kwargs):
    '''
    Ionize geometry via with supplied geometry and file containin list of ions.
    Parameters
    ----------
    geom : :obj:`~isicle.geometry.Geometry`
        Molecule representation.
    ion_method : str
        Alias for ionaztion method selection (explicit).
    **kwargs
        Keyword arguments to configure how mol objects can be saved.
        See :meth:`~isicle.adducts.ExplicitIonizationWrapper.finish`.
    Returns
    -------
    :obj:`~isicle.adducts.ExplicitIonizationWrapper`
        Dictionary of
    '''
    # Load ion file
    ion_list = load_ions(ion_path)

    # Parse ion file
    cations, anions, complex = parse_ions(ion_list)

    # Select ionization method
    iw = _ionize_method_selector(ion_method)

    # TODO make sure removeHs=False
    # Load in geometry of geom object
    iw.set_geomtery(geom)

    # Load specified ions by type
    iw.set_ions(cations, anions, complex)

    # TODO substructure check to validate specified ion support
    # iw.check_valid()

    # Generate adducts
    iw.generator()

    # Combine/ optionally write files
    res = iw.finish(**kwargs)

    return res


class ExplicitIonizationWrapper(IonizeInterface):
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

    def set_ions(self, cations: lst, anions: lst, complex: lst):
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

    def _add_ion(init_mol, base_atom_idx, ion_atomic_num):
        '''
        Adds specified ion (by atomic num) to specified base atom (by atom's index).
        Returns ionized RDKit mol object.
        '''
        import rdkit.Chem.rdchem as rdc
        from rdkit.Chem import AllChem
        # TODO only works for single atom ion

        # AddHs called here since not specified in load functions
        newmol = rdc.RWMol(init_mol)
        rdc.RWMol.AddAtom(newmol, rdc.Atom(ion_atomic_num))
        # atom idx starts at 0; GetNumAtoms returns last idx val+1
        atm_idx = init_mol.GetNumAtoms()
        rdc.RWMol.AddBond(newmol, base_atm_idx, atm_idx, Chem.BondType.SINGLE)
        Chem.SanitizeMol(newmol)
        # TODO insert try embed, need to catch any exceptions?
        AllChem.EmbedMolecule(newmol)
        AllChem.UFFOptimizeMolecule(newmol)
        return newmol

    def add_ion(init_mol, atom_dict, ion_atomic_num):
        '''
        Iterates through base atoms of RDKit mol object.
        Calls method to ionize and return new RDKit mol object.
        Returns dictionary of {base atom: ionized RDKit mol object}.
        '''

        mol_dict = {}
        for key in atom_dict.keys():
            newmol = self._add_ion(init_mol, key, ion_atomic_num)
            # Format {atom_base_idx:<newmol>}
            mol_dict[key] = newmol
        return mol_dict

    def positive_mode(self, init_mol, ion_atomic_num):
        '''
        Adds specified cation to RDKit mol object at every potential base atom.
        '''
        atom_dict = self.modify_atom_dict(init_mol, ion_atomic_num, mode='positive')
        mol_dict = self.add_ion(init_mol, atom_dict, ion_atomic_num)
        return mol_dict

    def _remove_ion(init_mol, base_atom, base_atom_idx, ion_atomic_num):
        '''
        Removes specified ion (by atomic num) from specified base atom (by atom's index).
        Returns deionized RDKit mol object.
        '''
        import rdkit.Chem.rdchem as rdc
        newmol = rdc.RWMol(init_mol)
        nbrs_dict = {nbr.GetAtomicNum(): nbr.GetIdx() for nbr in base_atom.GetNeighbors()}
        rdc.RWMol.RemoveAtom(newmol, nbrs_dict[ion_atomic_num])
        Chem.SanitizeMol(newmol)
        # TODO insert try embed, need to catch any exceptions?
        AllChem.EmbedMolecule(newmol)
        AllChem.UFFOptimizeMolecule(newmol)
        return newmol

    def remove_ion(init_mol, atom_dict, ion_atomic_num):
        '''
        Iterates through base atoms of RDKit mol object.
        Calls method to deionize and return new RDKit mol object.
        Returns dictionary of {base atom: deionized RDKit mol object}.
        '''

        mol_dict = {}
        for key, value in atom_dict.items():
            newmol = self._remove_ion(init_mol, key, value[3], ion_atomic_num)
            mol_dict[value[3]] = newmol
        return mol_dict

    def negative_mode(self, init_mol, ion_atomic_num):
        '''
        Removes specified anion from RDKit mol object at every potential base atom.
        '''
        atom_dict = self.modify_atom_dict(init_mol, mode='negative')
        mol_dict = self.remove_ion(init_mol, atom_dict, ion_atomic_num)
        return mol_dict

    def generator(self):
        '''
        Creates specified adducts at every atom site.
        '''
        pt = Chem.GetPeriodicTable()
        anions = self.anions
        cations = self.cations
        complex = self.complex
        init_mol = self.geom
        ion_dict = {}
        cation_dict = {}
        complex_dict = {}

        # TODO make sure removeHs=False in load functions, for now:
        init_mol = Chem.AddHs(init_mol)

        for ion in cations:
            sion = ion.strip('+')
            ion_atomic_num = pt.GetAtomicNumber(sion)
            mol_dict = self.positive_mode(init_mol, ion_atomic_num)
            cation_dict[ion] = mol_dict
            # cation_dict defined as {ion+:{base_atom_index: mol}}
        self.cations = cation_dict
        for ion in anions:
            sion = ion.strip('-')
            ion_atomic_num = pt.GetAtomicNumber(sion)
            mol_dict = self.negative_mode(init_mol, ion_atomic_num)
            anion_dict[ion] = mol_dict
            # anion_dict defined as {ion-:{base_atom_index: mol}}
        self.anions = anion_dict

        for ion in complex:
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


class CRESTIonizationWrapper(IonizeInterface):
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

    def set_ions(self, cations: lst, anions: lst, complex: lst):
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
        # consult https://github.com/grimme-lab/crest/issues/2 for supported ions by xtb CREST

    def generator(self, template=None, program='xtb', forcefield='gff', ewin=1, ion='H+', optlevel='Normal', dryrun=False):
        '''
        Call isicle.md.md for specified geom and every ion
        '''
        anions = self.anions
        cations = self.cations
        complex = self.complex
        ion_dict = {}
        cation_dict = {}
        complex_dict = {}
        # TODO check if md is returning mol object
        for x in cations:
            cation_dict[x] = md(self, program='xtb', template=template, tasks='protonate', forcefield=forcefield,
                                ewin=ewin, ion=x, optlevel=optlevel, dryrun=dryrun)
        self.cations = cation_dict
        for x in anions:
            anion_dict[x] = md(self, program='xtb', template=template, tasks='deprotonate', forcefield=forcefield,
                               ewin=ewin, ion=x, optlevel=optlevel, dryrun=dryrun)
        self.anions = anion_dict
        for x in complex:
            # TODO institute iterative approach to complex adduct formation
            # complex_dict[x] = None
            raise NotImplementedError
        self.complex = complex_dict

    def finish(self):
        '''
        Combine generated structures and optionally save

        Returns
        -------
        :obj:`~isicle.parse.NWChemResult`
            Parsed result data.
        '''
        ion_dict = {**self.cations, **self.anions, **self.complex}
        # ion dict format {ion<charge>:mol}
        return ion_dict
