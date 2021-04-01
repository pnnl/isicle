import abc


class FileParserInterface(metaclass=abc.ABCMeta):
    '''
    Abstract base class for file parser interface. All file parsers
    conform to this definition.
    '''

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'load')
                and callable(subclass.load)
                and hasattr(subclass, 'parse')
                and callable(subclass.parse)
                and hasattr(subclass, 'save')
                and callable(subclass.save)
                or NotImplemented)

    @abc.abstractmethod
    def load(self, path: str):
        '''Load in the data file'''
        raise NotImplementedError

    @abc.abstractmethod
    def parse(self):
        '''Extract relevant information from data'''
        raise NotImplementedError

    @abc.abstractmethod
    def save(self, path: str):
        '''Write parsed object to file'''
        raise NotImplementedError


class XYZGeometryInterface(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'dft_optimize')
                and callable(subclass.dft_optimize)
                and hasattr(subclass, 'md_optimize')
                and callable(subclass.md_optimize)
                and hasattr(subclass, 'get_natoms')
                and callable(subclass.get_natoms)
                and hasattr(subclass, 'get_atom_indices')
                and callable(subclass.get_atom_indices)
                and hasattr(subclass, 'get_global_properties')
                and callable(subclass.get_global_properties)
                and hasattr(subclass, '__copy__')
                and callable(subclass.__copy__)
                and hasattr(subclass, 'to_xyzblock')
                and callable(subclass.to_xyzblock)
                and hasattr(subclass, 'save_xyz')
                and callable(subclass.save_xyz)
                and hasattr(subclass, 'save_pickle')
                and callable(subclass.save_pickle)
                and hasattr(subclass, 'save')
                and callable(subclass.save)
                or NotImplemented)

    @abc.abstractmethod
    def dft_optimize(self):
        '''Optimize geometry using density function theory (DFT) methods.'''
        raise NotImplementedError

    @abc.abstractmethod
    def md_optimize(self):
        '''Optimize geometry using molecule dynamics methods (MD).'''
        raise NotImplementedError

    @abc.abstractmethod
    def get_natoms(self):
        '''Count number of atoms'''
        raise NotImplementedError

    @abc.abstractmethod
    def get_atom_indices(self):
        '''Extract indices of each atom from the internal geometry.'''
        raise NotImplementedError

    @abc.abstractmethod
    def get_global_properties(self):
        '''Return a copy of this object's global_properties dictionary'''
        raise NotImplementedError

    @abc.abstractmethod
    def __copy__(self):
        '''Return hard copy of this class instance.'''
        raise NotImplementedError

    @abc.abstractmethod
    def to_xyzblock(self):
        '''Get XYZ text for this structure.'''
        raise NotImplementedError

    @abc.abstractmethod
    def save_xyz(self):
        '''Save molecule as XYZ file'''
        raise NotImplementedError

    @abc.abstractmethod
    def save_pickle(self):
        '''Pickle this class instance.'''
        raise NotImplementedError

    @abc.abstractmethod
    def save(self, path: str):
        '''Write 3D molecule to file'''
        raise NotImplementedError


class GeometryInterface(XYZGeometryInterface):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'get_mol')
                and callable(subclass.get_mol)
                and hasattr(subclass, 'get_total_partial_charge')
                and callable(subclass.get_total_partial_charge)
                and hasattr(subclass, 'to_smiles')
                and callable(subclass.to_smiles)
                and hasattr(subclass, 'to_inchi')
                and callable(subclass.to_inchi)
                and hasattr(subclass, 'to_smarts')
                and callable(subclass.to_smarts)
                or NotImplemented)

    @abc.abstractmethod
    def to_mol(self, path: str):
        '''Returns RDKit Mol object for this Geometry.'''
        raise NotImplementedError

    @abc.abstractmethod
    def get_total_partial_charge(self):
        '''Determine total partial charge of molecule'''
        raise NotImplementedError

    @abc.abstractmethod
    def to_smiles(self, path: str):
        '''Return SMILES representation'''
        raise NotImplementedError

    @abc.abstractmethod
    def to_inchi(self, path: str):
        '''Return InChI representation'''
        raise NotImplementedError

    @abc.abstractmethod
    def to_smarts(self, path: str):
        '''Return SMARTS representation'''
        raise NotImplementedError


class AdductInterface(GeometryInterface):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'negative_mode')
                and callable(subclass.negative_mode)
                or NotImplemented)

    @abc.abstractmethod
    def negative_mode(self):
        '''Optimize geometry'''
        raise NotImplementedError


class WrapperInterface(metaclass=abc.ABCMeta):
    '''
    Abstract base class for wrapper interface. All QM
    wrappers conform to this definition.
    '''

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'load_geometry')
                and callable(subclass.load_geometry)
                and hasattr(subclass, 'set_geometry')
                and callable(subclass.set_geometry)
                and hasattr(subclass, 'configure')
                and callable(subclass.configure)
                and hasattr(subclass, 'save_config')
                and callable(subclass.save_config)
                and hasattr(subclass, 'run')
                and callable(subclass.run)
                and hasattr(subclass, 'finish')
                and callable(subclass.finish)
                or NotImplemented)

    @abc.abstractmethod
    def set_geometry(self):
        '''Load the input geometry file'''
        raise NotImplementedError

    @abc.abstractmethod
    def configure(self):
        '''Configure the run.'''
        raise NotImplementedError

    @abc.abstractmethod
    def save_config(self):
        '''Configure the run.'''
        raise NotImplementedError

    @abc.abstractmethod
    def run(self):
        '''Execute/submit the run.'''
        raise NotImplementedError

    @abc.abstractmethod
    def finish(self, path: str):
        '''Finalize, parse, return result object.'''
        raise NotImplementedError


class MDWrapperInterface(metaclass=abc.ABCMeta):
    '''
    Abstract base class for molecular dynamics wrapper interface. All QM
    wrappers conform to this definition.
    '''

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'load_geometry')
                and callable(subclass.load_geometry)
                and hasattr(subclass, 'set_geometry')
                and callable(subclass.set_geometry)
                and hasattr(subclass, 'job_type')
                and callable(subclass.job_type)
                and hasattr(subclass, 'run')
                and callable(subclass.run)
                and hasattr(subclass, 'finish')
                and callable(subclass.finish)
                or NotImplemented)

    @abc.abstractmethod
    def set_geometry(self):
        '''Load the input geometry file'''
        raise NotImplementedError

    @abc.abstractmethod
    def job_type(self):
        '''Make list of jobs.'''
        raise NotImplementedError

    @abc.abstractmethod
    def run(self):
        '''Execute/submit the run.'''
        raise NotImplementedError

    @abc.abstractmethod
    def finish(self, path: str):
        '''Finalize, parse, return result object.'''
        raise NotImplementedError
