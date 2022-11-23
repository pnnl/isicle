import os
import pickle
from io import StringIO

import isicle
import numpy as np
import pandas as pd
from isicle.interfaces import GeometryInterface, XYZGeometryInterface
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.SaltRemover import SaltRemover


class XYZGeometry(XYZGeometryInterface):
    '''
    Molecule information, specialized for XYZ files (bonding info not provided).
    It is not recommended to manipulate or retrieve
    attributes of this class without using class functions.

    Attributes
    ----------
    xyz: list of str
        Lines of XYZ block used to represent this compound's structure
    global_properties : dict
        Dictionary of properties calculated for this structure. Always at least
        contains the "from" key, which reports the last operation performed on
        this compound (e.g. load, dft_optimize). To generate compound
        properties, use get_* and *_optimize functions.
    history: list of str
        All steps performed on this compound since initial generation. For
        example, last step of history should always match "from".
    '''

    _defaults = ('history', 'xyz')
    _default_value = None

    def __init__(self, **kwargs):
        self.__dict__.update(dict.fromkeys(
            self._defaults, self._default_value))
        self.__dict__.update(kwargs)

        if self.history is None:
            self.history = []

    def _upgrade_to_Geometry(self, mol):
        '''
        Upgrade this class to the Geometry class using a provided mol object.
        Hard copies everything but the xyz structure.

        Parameters
        ----------
        mol : RDKit Mol object
            Structure to use in new Geometry object

        Returns
        -------
        Geometry
            Instance containing the new mol representation.

        '''

        # Create dict to load in
        d = self.__dict__.copy()
        d['history'] = self.get_history()
        d['mol'] = mol
        d.pop('xyz')

        return Geometry(**d)

    def _update_structure(self, inplace, mol=None, xyz=None, xyz_filename=None,
                          event=None):
        '''
        Return object with given structure. If inplace and a xyz structure is
        provided, manipulates the current instance. Otherwise, a new
        XYZGeometry or Geometry object is created with hard copied attributes.
        Only one of mol, xyz, or xyz_filename is to be provided. If multiple
        are provided, they are prioritized in the following order: (1) mol,
        (2) xyz_filename, (3) xyz. If none of these are provided, an error is
        raised.

        Parameters
        ----------
        inplace : boolean
            If true, update this instance with the new structure. Otherwise,
            create a new XYZGeometry or Geometry instance and populate it with
            the structure.
        mol : RDKit Mol object (optional)
            Structure with 3D abilities.
        xyz_filename: str (optional)
            Path to file to load in xyz from
        xyz: list of str (optional)
            Lines from xyz block to use as structure representation.
        event: str (optional)
            Operation that led to this structure change, for history-tracking
            purposes within the instance returned.

        Returns
        -------
        XYZGeometry or Geometry
            Updated structure instance. Geometry is returned if mol was
            provided, otherwise XYZGeometry is returned.

        '''

        if mol is not None:
            # Upgrade to the Geometry class
            geom = self._upgrade_to_Geometry(mol)

        else:
            if inplace:  # Modify this object
                geom = self
            else:  # Make a new object
                geom = self.__copy__()

            # Ensure exactly one of xyz or a filename is provided
            # Add as structure of the new XYZGeometry
            if xyz_filename is not None:
                geom.xyz = load_xyz(xyz_filename)
            elif xyz is not None:
                geom.xyz = xyz

        # Add event that led to this change in structure
        # If event is None, nothing will happen
        geom._update_history(event)

        return geom

    def _update_history(self, event):
        '''Add event to this instance's history and update 'from' in the
        global properties dictionary. If none, nothing is updated. Returns
        a copy of the full history.'''
        if event is not None:
            self.history.append(event)
            self.__dict__['from'] = event
        return self.get_history()

    def get_history(self):
        '''Returns a copy of this object's history (events that led to its
        current state).'''
        return self.history[:]

    def get_basename(self):
        '''Returns a copy of this object's basename (original filename).'''
        return self.basename

    def add___dict__(self, d, override=False):
        '''Accepts a dictionary of values and adds any non-conflicting
        information to the attribute dictionary'''

        if not override:
            # Remove keys that are already present in attribute dictionary
            [d.pop(k, None) for k in self.__dict__.keys()]

        # Add to attributes
        self.__dict__.update(d)
        return

    def dft(self, program='NWChem', template=None, **kwargs):
        '''
        Optimize geometry from XYZ, using stated functional and basis set.
        Additional inputs can be grid size, optimization criteria level,
        '''
        geom = isicle.qm.dft(self.__copy__(), program=program,
                             template=template, **kwargs)

        return geom

    def md(self, program='xtb', **kwargs):
        '''
        Optimize geometry or generate conformers or adducts from XYZ using stated forcefield.
        Additional inputs can be energy window, optimization criteria level, charge, or ion.
        '''
        geom = isicle.md.md(self.__copy__(), program=program, **kwargs)

        return geom

    def ionize(self, ion_path=None, ion_list=None, save=False, path=None, **kwargs):
        '''
        Calls xtb CREST to ionize xyz geometry, using specified list of ions and method of ionization.
        WARNING: must run isicle.geometry.XYZGeometry.set_formal_charge prior to this.

        Parameters
        ----------
        ion_path : str
            Filepath to text file containing ions with charge (eg. `H+`) to be considered
            Either ion_path or ion_list must be specified
        ion_list : list
            List of strings of adducts to be considered.
            Must be specifed in syntax `Atom+` or `Atom-`, eg. `H+`, `Na+`, `H-Na+`
            Either ion_path or ion_list must be specified
        save : bool
            Saves wrapper object to .pkl in specified path directory
        path : str (optional)
            Directory to write output files
        **kwargs :
            see :meth: `~isicle.adducts.CRESTIonizationWrapper.submit`
            for more options

        Returns
        -------
        Dictionary of adducts, `{<IonCharge>:[<geomObjects>]}`
        '''
        if self.__dict__.get('charge') is None:
            raise ValueError(
                'Must first run isicle.geometry.XYZGeometry.set_formal_charge for an xyz structure')
        iw = isicle.adducts.ionize("crest").run(
            self.__copy__(), ion_path=ion_path, ion_list=ion_list, **kwargs)

        if save == True and path is not None:
            iw.save(path)

        return iw

    def set_formal_charge(self, charge):
        '''Specify and set the known formal charge of the xyz molecule to __dict__'''
        try:
            charge = int(charge)
        except:
            raise TypeError('Charge specifed must be str or int')
        self.__dict__.update(charge=charge)
        return self.charge

    def get_natoms(self):
        '''Calculate total number of atoms.'''
        self.__dict__['natoms'] = int(self.xyz[0].strip())
        return self.__dict__['natoms']

    def get_atom_indices(self, atoms=['C', 'H']):
        '''
        Extract indices of each atom from the internal geometry.

        Parameters
        ----------
        atoms : list of str
            Atom types of interest.

        Returns
        -------
        list of int
            Atom indices.
        '''
        idx = []

        for i in range(2, len(self.xyz)):
            atom = self.xyz[i].split(' ')[0]
            if atom in atoms:
                idx.append(i - 2)
        return idx

    def get___dict__(self):
        '''Return a copy of this object's attribute dictionary'''
        return self.__dict__.copy()

    def __copy__(self):
        '''Return hard copy of this class instance.'''
        # TODO: manage what should be passed, rather than all?
        d = self.__dict__.copy()
        return type(self)(**d)

    def to_xyzblock(self):
        '''Get XYZ text for this structure.'''
        return '\n'.join(self.xyz)

    def save_xyz(self, path):
        '''Save molecule as XYZ file'''
        with open(path, 'w') as f:
            f.write(self.to_xyzblock())
        return 'Success'

    def save_pickle(self, path):
        '''Pickle this class instance.'''
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        return 'Success'

    def save_mfj(self, path):
        # Extract properties
        props = self.get___dict__()

        # Check for charges in global properties
        if ('energy' not in props) or ('charge' not in props):
            raise KeyError('DFT energy optimization required. '
                           'See isicle.qm.dft.')

        # Get XYZ coordinates
        xyz = pd.read_csv(StringIO(self.to_xyzblock()),
                          skiprows=2, header=None,
                          delim_whitespace=True,
                          names=['Atom', 'x', 'y', 'z'])

        # Extract and append charges
        xyz['Charge'] = props['charge']

        # Load masses and merge
        masses = isicle.utils.atomic_masses()[['Symbol', 'Mass']]
        mfj = pd.merge(xyz, masses, left_on='Atom', right_on='Symbol')

        # Rename columns
        mfj = mfj[['x', 'y', 'z', 'Mass', 'Charge']].astype(float)

        # Write to file
        with open(path, 'w') as f:
            f.write(os.path.splitext(os.path.basename(path))[0] + '\n')
            f.write('1\n')
            f.write(str(len(mfj.index)) + '\n')
            f.write('ang\n')
            f.write('calc\n')
            f.write('1.000\n')

            for row in mfj.values:
                f.write('\t'.join([str(x) for x in row]) + '\n')

    def save(self, path, fmt=None):
        '''
        Save molecule

        Parameters
        ----------
        path : str
            Path to save file.
        fmt : str (optional)
            Format to save this molecule in. If None, determined from given
            path's extension. If .pkl. pickles this full object.
            Default: None.
            Supported formats: .xyz, .pkl.

        Returns
        -------
        str
            Status of save.

        '''

        if fmt is None:
            # Decide format based on path
            fmt = os.path.splitext(path)[-1]
        fmt = fmt.lower()

        if 'xyz' in fmt:
            return self.save_xyz(path)

        if 'pkl' in fmt:
            return self.save_pickle(path)

        if 'mfj' in fmt:
            return self.save_mfj(path)

        raise TypeError(
            'Input format {} not supported for {}.'.format(fmt, self.__class__))


class Geometry(XYZGeometry, GeometryInterface):
    '''
    Molecule information, specialized for 3D representations.
    It is not recommended to manipulate or retrieve
    attributes of this class without using class functions.

    Attributes
    ---------- =
    mol : RDKit Mol object
        Current structure.
    __dict__ : dict
        Dictionary of properties calculated for this structure.
    history: list of str
        All steps performed on this compound since initial generation. For
        example, last step of history should always match "from".
    '''

    _defaults = ('history', 'mol')
    _default_value = None

    def __init__(self, **kwargs):
        self.__dict__.update(dict.fromkeys(
            self._defaults, self._default_value))
        self.__dict__.update(kwargs)

        if self.history is None:
            self.history = []

    def _downgrade_to_XYZGeometry(self, xyz):
        '''
        Downgrade this class to the XYZGeometry class using a provided xyz.
        Hard copies everything but the mol structure.

        Parameters
        ----------
        xyz : list of str
            Lines of xyz block to use for new object.

        Returns
        -------
        XYZGeometry
            Instance containing the new xyz representation.

        '''

        # Create dict to load in
        d = self.__dict__.copy()
        d['xyz'] = xyz
        d.pop('mol')

        return XYZGeometry(**d)

    def _update_structure(self, inplace, mol=None, xyz=None, xyz_filename=None,
                          event=None):
        '''
        Return object with given structure. If inplace and a mol structure is
        provided, manipulates the current instance. Otherwise, a new
        XYZGeometry or Geometry object is created with hard copied attributes.
        Only one of mol, xyz, or xyz_filename is to be provided. If multiple
        are provided, they are prioritized in the following order: (1) mol,
        (2) xyz_filename, (3) xyz. If none of these are provided, an error is
        raised.

        Parameters
        ----------
        inplace : boolean
            If true, update this instance with the new structure. Otherwise,
            create a new XYZGeometry or Geometry instance and populate it with
            the structure.
        mol : RDKit Mol object (optional)
            Structure with 3D abilities.
        xyz_filename: str (optional)
            Path to file to load in xyz from
        xyz: list of str (optional)
            Lines from xyz block to use as structure representation.
        event: str (optional)
            Operation that led to this structure change, for history-tracking
            purposes within the instance returned.

        Returns
        -------
        XYZGeometry or Geometry
            Updated structure instance. Geometry is returned if mol was
            provided, otherwise XYZGeometry is returned.

        '''

        if mol is None:

            if xyz_filename is not None:
                # Load in file first
                xyz = load_xyz(xyz_filename)

            if xyz is not None:
                # Downgrade to XYZGeometry class
                geom = self._downgrade_to_XYZGeometry(xyz)

        else:
            if inplace:  # Modify this object
                geom = self
            else:  # Make a new object
                geom = self.__copy__()
            geom.mol = mol

        # Add event that led to this change in structure
        # If event is None, nothing will happen
        geom._update_history(event)

        return geom

    def initial_optimize(self, embed=False, forcefield='UFF', ff_iter=200, inplace=False):
        '''
        Params
        ------
        embed: bool
            Try embedding molecule with eigenvalues of distance matrix
            Except failure seeds with random coordinates
        forcefield: str
            Specify "UFF" for Universal Force Field optimization by RDKit
            Specify "MMFF" or "MMFF94" for Merck Molecular Force Field 94
            Specify "MMFF94s" for the MMFF94 s variant
        inplace: bool
            If geom.mol is modified inplace
        Returns
        -------
        Geometry
            Updated structure instance
        '''
        def _forcefield_selector(forcefield, mw):
            '''
            Checks if user specified forcefield is compatible with supplied mol
            '''
            forcefield = forcefield.lower()
            if forcefield == 'uff':
                if Chem.rdForceFieldHelpers.UFFHasAllMoleculeParams(mw) is True:
                    return Chem.AllChem.UFFOptimizeMolecule
                else:
                    raise ValueError(
                        'UFF is not available for all atoms in molecule.')
            elif forcefield in ['mmff', 'mmff94', 'mmff94s']:
                if Chem.rdForceFieldHelpers.MMFFHasAllMoleculeParams(mw) is True:
                    return Chem.rdForceFieldHelpers.MMFFOptimizeMolecule
                else:
                    raise ValueError(
                        'specified MMFF is not available for all atoms in molecule.')
            else:
                raise ValueError(
                    'RDKit only supports UFF, MMFF, MMFF94, MMFF94s as forcefields.')

        mol = self.to_mol()
        if embed == True:
            try:
                Chem.AllChem.EmbedMolecule(mol)
            except:
                Chem.AllChem.EmbedMolecule(mol, useRandomCoords=True)
        if forcefield is not None:
            if 'mmff94s' in forcefield.lower():
                _forcefield_selector(forcefield, mol)(
                    mol, mmffVariant=forcefield, maxIters=ff_iter)
            else:
                _forcefield_selector(forcefield, mol)(mol, maxIters=ff_iter)

        geom = self._update_structure(inplace, mol=mol)
        geom._update_history('initial_optimize')

        return geom

    def desalt(self, salts=None, inplace=False):
        '''
        Desalts RDKit mol object using Chem.SaltRemover module.

        Parameters
        ----------
        salts : str (optional)
            Salt type to remove. Ex: 'Cl', 'Br', '[Na+]'. Default: None.
        inplace : boolean (optional)
            If true, update this instance with the new structure. Otherwise,
            create a new Geometry instance and populate it with the structure.
            Default: False.

        Returns
        -------
        Geometry
            Molecule with given salt(S) removed.

        '''

        # TODO: should this raise an error instead?
        # If no salts given, skip desalting
        if salts is None:
            return self._update_structure(inplace, mol=self.to_mol())

        remover = SaltRemover(defnFormat='smiles', defnData=salts)
        # defnData="[Cl,Br,Na]" *sample definition of salts to be removed*
        # add iterator for salts listed in config?
        # set potential salts to be removed in a config file

        mol, deleted = remover.StripMolWithDeleted(self.to_mol())
        # using StripMolWithDeleted instead of StripMol
        # add functionality to track removed salts
        # atomno = res.GetNumAtoms
        # if relevant to future use, returns atom count post desalting

        geom = self._update_structure(inplace, mol=mol)
        geom._update_history('desalt')
        # TODO: add any properties from this operation to global_props?

        return geom

    def neutralize(self, inplace=False):
        '''
        Neutralizes RDKit mol object using neutralization reactions.

        Parameters
        ----------
        inplace : boolean (optional)
            If true, update this instance with the new structure. Otherwise,
            create a new Geometry instance and populate it with the structure.
            Default: False.

        Returns
        -------
        Geometry
            Neutralized form of the molecule.

        '''

        def _initialize_neutralisation_reactions():
            patts = (
                # Imidazoles
                ('[n+;H]', 'n'),
                # Amines
                ('[N+;!H0]', 'N'),
                # Carboxylic acids and alcohols
                ('[$([O-]);!$([O-][#7])]', 'O'),
                # Thiols
                ('[S-;X1]', 'S'),
                # Sulfonamides
                ('[$([N-;X2]S(=O)=O)]', 'N'),
                # Enamines
                ('[$([N-;X2][C,N]=C)]', 'N'),
                # Tetrazoles
                ('[n-]', '[nH]'),
                # Sulfoxides
                ('[$([S-]=O)]', 'S'),
                # Amides
                ('[$([N-]C=O)]', 'N'),
            )
            return [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y)) for x, y in patts]

        reactions = _initialize_neutralisation_reactions()

        mol = self.to_mol()
        replaced = False
        for i, (reactant, product) in enumerate(reactions):
            while mol.HasSubstructMatch(reactant):
                replaced = True
                rms = Chem.AllChem.ReplaceSubstructs(mol, reactant, product)
                mol = rms[0]

        geom = self._update_structure(inplace, mol=mol)
        geom._update_history('neutralize')
        # TODO: add any properties from this operation to global_props?

        return geom

    def tautomerize(self, return_all=False, inplace=False):
        '''
        Generate tautomers according to RDKit TautomerEnumerator() method.

        Parameters
        ----------
        return_all : boolean (optional)
            If true, return all tautomers generated. Otherwise, only return
            the most common. Default=False
        inplace : boolean (optional)
            If true, update this instance with the new structure. Otherwise,
            create a new Geometry instance and populate it with the structure.
            Default: False.

        Returns
        -------
        Geometry or list(Geometry)
            Tautomer(s) of starting structure.

        '''

        # source: https://rdkit.blogspot.com/2020/01/trying-out-new-tautomer.html
        # Discuss noted double bond changes
        enumerator = rdMolStandardize.TautomerEnumerator()

        mol = self.to_mol()
        res = [mol]
        tauts = enumerator.Enumerate(mol)
        smis = [Chem.MolToSmiles(x) for x in tauts]
        s_smis = sorted((x, y)
                        for x, y in zip(smis, tauts) if x != self.to_smiles())
        res += [y for x, y in s_smis]

        # Ensure res is a list of mol objects
        if return_all:
            new_geoms = []
            for r in res:
                geom = self._update_structure(
                    False, mol=r, event='tautomerize')
                new_geoms.append(geom)
            return new_geoms

        geom = self._update_structure(inplace, mol=res[0])
        geom._update_history('tautomerize')
        # TODO: add any properties from this operation to global_props?

        return geom

    def kekulize(self, inplace=False):
        '''
        `Kekulizes the molecule if possible. Otherwise the molecule is not modified`
        This is recommended for aromatic molecules undergoing explicit ionization.
        Aromatic bonds are explicitly defined and aromatic flags are cleared.
        '''
        mol = Chem.KekulizeIfPossible(self.to_mol(), clearAromaticFlags=True)
        geom = self._update_structure(inplace, mol=mol, event='kekulize')
        geom._update_history('kekulize')
        return geom

    def ionize(self, ion_path=None, ion_list=None, method='explicit', save=False, path=None, **kwargs):
        '''
        Ionize geometry, using specified list of ions and method of ionization.

        Parameters
        ----------
        ion_path : str
            Filepath to text file containing ions with charge (eg. `H+`) to be considered
            Either ion_path or ion_list must be specified
        ion_list : list
            List of strings of adducts to be considered.
            Must be specifed in syntax `Atom+` or `Atom-`, eg. `H+`, `Na+`, `H-Na+`
            Either ion_path or ion_list must be specified
        method : str
            Method of ionization to be used, 'explicit' or 'crest' is accepted
        save : bool
            Saves wrapper object to .pkl in specified path directory
        path : str (optional)
            Directory to write output files
        ensembl : bool (optional)
            Returns instead a list of adduct geometries
        **kwargs:
            see :meth: `~isicle.adducts.ExplicitIonizationWrapper.submit` or
                       `~isicle.adducts.CRESTIonizationWrapper.submit`
            for more options

        Returns
        -------
        Dictionary of adducts, `{<IonCharge>:[<geomObjects>]}`
        '''
        iw = isicle.adducts.ionize(method).run(
            self.__copy__(), ion_path=ion_path, ion_list=ion_list, **kwargs)

        if save == True and path is not None:
            iw.save(path)

        return iw

    def get_natoms(self):
        '''Calculate total number of atoms.'''
        natoms = Chem.Mol.GetNumAtoms(self.to_mol())
        self.__dict__['natoms'] = natoms
        return self.__dict__['natoms']

    def get_atom_indices(self, atoms=['C', 'H'], lookup={'C': 6, 'H': 1, 'N': 7,
                                                         'O': 8, 'F': 9, 'P': 15}):
        '''
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
        '''

        atoms = [lookup[x] for x in atoms]
        idx = []
        for a in self.mol.GetAtoms():
            if a.GetAtomicNum() in atoms:
                idx.append(a.GetIdx())

        return idx

    def get_total_partial_charge(self):
        '''Sum the partial charge across all atoms.'''
        mol = self.to_mol()
        Chem.AllChem.ComputeGasteigerCharges(mol)
        contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge')
                    for i in range(mol.GetNumAtoms())]
        return np.nansum(contribs)

    def get_formal_charge(self):
        '''Calculate formal charge of the molecule.'''
        mol = self.to_mol()
        return Chem.rdmolops.GetFormalCharge(mol)

    def set_formal_charge(self):
        '''Set formal charge of the molecule to __dict__'''
        self.__dict__.update(charge=self.get_formal_charge())
        return self.charge

    def __copy__(self):
        '''Return hard copy of this class instance.'''
        # TODO: manage what should be passed, rather than all?
        d = self.__dict__.copy()
        d['history'] = self.get_history()
        d['mol'] = self.to_mol()
        return type(self)(**d)

    def to_mol(self):
        '''
        Returns RDKit Mol object for this Geometry.

        Returns
        -------
        RDKit Mol object
            Copy of structure

        '''

        return self.mol.__copy__()

    def to_smiles(self):
        '''Get SMILES for this structure.'''
        return Chem.MolToSmiles(self.to_mol())

    def to_inchi(self):
        '''Get InChI for this structure.'''
        return Chem.MolToInchi(self.to_mol())

    def to_smarts(self):
        '''Get SMARTS for this structure.'''
        return Chem.MolToSmarts(self.to_mol())

    def to_xyzblock(self):
        '''Get XYZ text for this structure.'''
        return Chem.MolToXYZBlock(self.to_mol())

    def to_pdbblock(self):
        '''Get PDB text for this structure'''
        return Chem.MolToPDBBlock(self.to_mol())

    def to_molblock(self):
        '''Get PDB text for this structure'''
        return Chem.MolToMolBlock(self.to_mol())

    def save_smiles(self, path: str):
        '''Save this structure's SMILES to file.'''
        with open(path, 'w') as f:
            f.write(self.to_smiles())
        return 'Success'

    def save_inchi(self, path: str):
        '''Save this structure's InChI to file.'''
        with open(path, 'w') as f:
            f.write(self.to_inchi())
        return 'Success'

    def save_smarts(self, path: str):
        '''Save this structure's SMARTS to file.'''
        with open(path, 'w') as f:
            f.write(self.to_smarts())
        return 'Success'

    def save_xyz(self, path: str):
        '''Save XYZ file for this structure.'''
        return Chem.MolToXYZFile(self.to_mol(), path)

    def save_mol(self, path):
        '''Save Mol file for this structure.'''
        return Chem.MolToMolFile(self.to_mol(), path)

    def save_pdb(self, path: str):
        '''Save PDB file for this structure.'''
        return Chem.MolToPDBFile(self.mol, path)

    # def save_mfj(self, path):
    #     self._downgrade_to_XYZGeometry(self.to_xyzblock()).save_mfj(path)

    def save(self, path, fmt=None):
        '''
        Save molecule

        Parameters
        ----------
        path : str
            Path to save file.
        fmt : str (optional)
            Format to save this molecule in. If None, determined from given
            path's extension. If .pkl. pickles this full object.
            Default: None.
            Supported formats: .smi (SMILES), .inchi (InChI), .mol, .xyz,
            .pdb, .pkl.

        Returns
        -------
        str
            Status of save.

        '''

        if fmt is None:
            # Decide format based on path
            fmt = os.path.splitext(path)[-1]
        fmt = fmt.lower()

        if 'smi' in fmt:
            return self.save_smiles(path)

        if 'inchi' in fmt:
            return self.save_inchi(path)

        if 'smarts' in fmt:
            return self.save_smarts(path)

        if 'xyz' in fmt:
            return self.save_xyz(path)

        if 'mol' in fmt:  # diff for mol2?
            return self.save_mol(path)

        if 'pkl' in fmt:
            return self.save_pickle(path)

        if 'pdb' in fmt:
            return self.save_pdb(path)

        if 'mfj' in fmt:
            return self.save_mfj(path)

        # TODO: enable Compute2DCoords, https://www.rdkit.org/docs/source/rdkit.Chem.rdDepictor.html

        raise TypeError(
            'Input format {} not supported for {}.'.format(fmt, self.__class__))


class MDOptimizedGeometry(Geometry):
    '''
    Builds off of the 3D representation, with additional methods specific to a
    representation with MD optimized 3D coordinates. Any methods that would
    result in a more defined representation (e.g. DFT optimized) should yield
    the appropriate subclass.
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()


class DFTOptimizedGeometry(Geometry):
    '''
    Builds off of the 3D representation, with additional methods specific to a
    representation with DFT optimized 3D coordinates.
    '''

    def __init__(self):
        # Initialize the base class
        super().__init__()
