from isicle.interfaces import AdductInterface
from isicle import geometry as geom
import os
import rdkit.Chem.rdchem as rdc
from rdkit import Chem


# TODO institute adduct calls in geometry


def load_ions(self, path: str):
    '''
    Read adduct ions from file.

    Parameters
    ----------
    path : str
        Path to ions text file (.txt, etc.), each ion should be a new line.

    Returns
    -------
    list
        List of ions from given text file.
    '''
    #self.path = path
    contents = geom._load_text(path)
    self.contents = [line.rstrip() for line in contents]
    return self.contents


def _select_method(method):
    method_map = {'explicit': ExplicitIonizationMachine, 'crest': CRESTIonizationWrapper}

    if method.lower() in method_map.keys():
        return method_map[method.lower()]()
    else:
        raise ValueError('{} not a supported ionization method.'.format(method))


def ionize(self, path: str, ionpath: str, method: str):
    '''
    '''
    # Load geometry
    geom = geom.load(path)

    # Instiate Adduct instance
    adds = Adduct()

    # Populate instance information
    adds.path = geom.path

    # Load ion file
    adds.contents = load_ions(ionpath)

    # Parse ions
    adds.parse_ions()

    # Select ionization method
    adds.method = _select_method(method)

    # Generate adducts and select method
    adds.call_method()

    # Write returned structures to file
    # adds.save()


class Adduct():

    def __init__(self, mol=None, anions=None, cations=None, complex=None, path=None, contents=None, method=None):
        self.mol = mol
        self.anions = anions
        self.cations = cations
        self.complex = complex
        self.path = path
        self.contents = contents
        self.method = method

    def parse_ions(self, anions=None, cations=None, complex=None):
        '''
        Parse and categorize list of ion by ion type.

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
        if anions is None:
            anions = []
        if cations is None:
            cations = []
        if complex is None:
            complex = []

        ions = self.contents
        # TODO make parser more robust for input list variation
        for x in ions:
            if ('+' and '-') in x:
                complex.append(x)
            elif ('+') in x:
                cations.append(x)
            elif ('-') in x:
                anions.append(x)
            else:
                raise ValueError('Unrecognized ion specified.')
                # TODO insert more descriptive error

        self.anions = anions
        self.cations = cations
        self.complex = complex

    def call_method(self):
        '''
        Call specified ionization class.
        '''
        method = self.method

        return method()

    def save(self):
        '''
        Writes output molecules to file.
        '''
        # TODO make options to save / write from finish methods in specified ion method


class ExplicitIonizationMachine(IonizeInterface):
    def __init__(self, mol_dict=None):
        mol_dict = None

    def generator(self):
        '''
        Creates specified adducts at every atom site.
        '''
        pt = Chem.GetPeriodicTable()
        anions = self.anions
        cations = self.cations
        complex = self.complex
        init_mol = self.mol
        anion_dict = {}
        cation_dict = {}
        complex_dict = {}
        # self.check_valid()
        # TODO substructure check to validate specified ion support

        for ion in anions:
            ion_atomic_num = pt.GetAtomicNumber(ion)
            mol_dict = self.negative_mode(init_mol, ion_atomic_num)
            anion_dict[ion] = mol_dict
            # anion_dict defined as {ion:{base_atom_index: mol}}

        for ion in cations:
            ion_atomic_num = pt.GetAtomicNumber(ion)
            mol_dict = self.positive_mode(init_mol, ion_atomic_num)
            cation_dict[ion] = mol_dict
            # cation_dict defined as {ion:{base_atom_index: mol}}

        for ion in complex:
            # TODO iterative method for adduct generation
            # Call both positve_mode and negative_mode
            raise NotImplementedError

    def check_valid(self):
        '''
        Performs substructure search and checks if proposed adduct can be formed.
        '''
        # consider calc_new_mz from
        # https://github.com/pnnl/mame/blob/24f9fcb19639d5a1c2ca0f788c55d6a2efe18ca6/mame/mame_utils.py
        # consider add_formala, remove_formula from
        # https://github.com/pnnl/mame/blob/24f9fcb19639d5a1c2ca0f788c55d6a2efe18ca6/mame/formula_module.py
        mol.HasSubstructMatch((ion.strip('+')).strip('-'))
        # complex structure checking, spectre?
        # consider how substructure would break off, not all atoms may be from same region

        # TODO (1) check if ions are in supported list of ions
        # TODO (2) modify self.anions, self.cations, self.complex to conform to supported ions

    def modify_atom_dict(self, init_mol, mode=None):
        all_atom_dict = {atom: [atom.GetSymbol(), atom.GetTotalValence() -
                                atom.GetDegree(), atom.GetTotalNumHs(includeNeighbors=True), atom.GetIdx()] for atom in init_mol.GetAtoms()}
        if mode == 'positive':
            atom_dict = {}
            for key, value in all_atom_dict.items():
                if value[0] in self.cations:
                    # TODO add check for more ion types
                    # eg. if value[0] in CationList
                    continue
                elif value[1] == 0 and value[2] == 0:
                    # removes atoms that cannot accept a bond
                    continue
                atom_dict[key] = value
            return atom_dict

        elif mode == 'negative':
            atom_dict = {}
            for key, value in all_atom_dict.items():
                if value[0] in self.anions:
                    # TODO add check for more ion types
                    # eg. if value[0] in AnionList
                    continue
                atom_dict[key] = value
            return atom_dict

        if mode == None:
            raise NotImplementedError

    def positive_mode(self, init_mol, ion_atomic_num):
        '''
        Adds specified cation to RDKit mol object at every potential base atom.
        '''
        atom_dict = self.modify_atom_dict(init_mol, mode='positive')
        mol_dict = self.add_ion(init_mol, atom_dict, ion_atomic_num)
        return mol_dict

    def _add_ion(init_mol, base_atom_idx, ion_atomic_num):
        '''
        Adds specified ion (by atomic num) to specified base atom (by atom's index).
        Returns ionized RDKit mol object.
        '''
        newmol = rdc.RWMol(init_mol)
        atm = rdc.Atom(ion_atomic_num)
        atm_idx = mol.GetNumAtoms()
        rdc.RWMol.AddAtom(newmol, atm)
        rdc.RWMol.AddBond(newmol, base_atm_idx, atm_idx, Chem.BondType.SINGLE)
        return newmol

    def add_ion(init_mol, atom_dict, ion_atomic_num):
        '''
        Iterates through base atoms of RDKit mol object.
        Calls method to ionize and return new RDKit mol object.
        Returns dictionary of {base atom: ionized RDKit mol object}.
        '''

        mol_dict = {}
        for key, value in atom_dict.items():
            newmol = self._add_ion(init_mol, key, ion_atomic_num)
            mol_dict[key] = newmol
        return mol_dict

    def negative_mode(self, init_mol, ion_atomic_num):
        '''
        Removes specified anion from RDKit mol object at every potential base atom.
        '''
        atom_dict = self.modify_atom_dict(init_mol, mode='negative')
        mol_dict = self.remove_ion(init_mol, atom_dict, ion_atomic_num)
        return mol_dict

    def _remove_ion(init_mol, base_atom, base_atom_idx, ion_atomic_num):
        '''
        Removes specified ion (by atomic num) from specified base atom (by atom's index).
        Returns deionized RDKit mol object.
        '''
        newmol = rdc.RWMol(init_mol)
        nbrs_dict = {pt.GetAtomicNumber(nbr): nbr.GetIdx() for nbr base_atom.GetNeighbors()}
        if ion_atomic_num in nbrs_dict.keys():
            rdc.RWMol.RemoveAtom(newmol, nbrs_dict[ion_atomic_num])
            bond_idx = newmol.GetBondIdx(base_atom_idx, nbrs_dict[ion_atomic_num])
            rdc.RWMol.RemoveBond(newmol, bond_idx)
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

    def finish(self):
        '''
        Writes finalized output mol blocks in dictionary to output files
        '''
        # TODO formalize this method


class CRESTIonizationWrapper(IonizeInterface):
    def __init__(self, path=path, cmd=None, ions=None, alt=False, custom=False):
        super().__init__()
        self.path = path
        self.cmd = cmd
        self.anion = anions

    def generation(self):
        '''
        Alternative adduct generation method.
        Utilizes xtb CREST software from Grimme group.
        Documentation:
        '''
        self.positive_mode()
        self.negative_mode()
        # TODO institute iterative approach to complex adduct formation

    def check_valid(self, ion):
        '''
        Performs substructure search and checks if supported by CREST documentation
        '''
        mol.HasSubstructMatch((ion.strip('+')).strip('-'))
        # TODO institute ion list check against list of supported CREST ions
        # consult https://github.com/grimme-lab/crest/issues/2 for supported ions by xtb CREST

    def positive_mode(self, path, cation=None, alt=False, custom=False):
        '''
        Generate xtb CREST protonate tool command
        https://xtb-docs.readthedocs.io/en/latest/crestxmpl.html#protonation-site-screening

        Default:
        [H+]: (ex.) crest 'struc.xyz' -protonate

        Alternate:
        (ex.) crest 'struc.xyz' -protonate -swel STR
        STR must be string of element and charge.
        [Na+]: (ex.) crest 'struc.xyz' -protonate -swel 'Na+'
        '''
        # consider -twel for tautomeric forms of output adduct
        path = self.path
        if self.alt:
            self.cmd = 'crest {path} -protonate -swel {cation}'.format(path=path, cation=cation)
        elif self.custom:
            # TODO parse custom config file
            raise NotImplementedError
        else:
            self.cmd = 'crest {path} -protonate'.format(path=path)
        return self

    def negative_mode(self, path, anion=None, alt=False, custom=False):
        '''
        Generate xtb CREST deprotonate tool command

        https://xtb-docs.readthedocs.io/en/latest/crestxmpl.html#deprotonation-site-screening

        Default:
        [H-]: (ex.) crest 'struc.xyz' -deprotonate
        '''
        # consider -twel for tautomeric forms of output adduct
        path = self.path
        if self.alt:
            self.cmd = 'crest {path} -deprotonate -swel {anion}'.format(path=path, anion=anion)
        elif self.custom:
            # TODO parse custom config file
            raise NotImplementedError
        else:
            self.cmd = 'crest {path} -deprotonate'.format(path=path)
        return self

    def run(self):
        '''
        Run xtb CREST protonate/deprotonate tool command
        '''
        os.system(self.cmd)
        # instititute temp dir

    def finish(self):
        '''

        '''
        # compile structures
