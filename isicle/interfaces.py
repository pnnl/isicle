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


class GeometryInterface(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'optimize')
                and callable(subclass.optimze)
                and hasattr(subclass, 'save')
                and callable(subclass.save)
                and hasattr(subclass, 'save')
                and callable(subclass.save)
                and hasattr(subclass, 'to_smiles')
                and callable(subclass.to_smiles)
                and hasattr(subclass, 'to_inchi')
                and callable(subclass.to_inchi)
                and hasattr(subclass, 'to_smarts')
                and callable(subclass.to_smarts)
                and hasattr(subclass, 'natoms')
                and callable(subclass.natoms)
                and hasattr(subclass, 'total_partial_charge')
                and callable(subclass.total_partial_charge)
                or NotImplemented)

    @abc.abstractmethod
    def optimize(self):
        '''Optimize geometry'''
        raise NotImplementedError

    @abc.abstractmethod
    def save(self, path: str):
        '''Write 3D molecule to file'''
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

    @abc.abstractmethod
    def natoms(self):
        '''Count number of atoms'''
        raise NotImplementedError

    @abc.abstractmethod
    def total_partial_charge(self):
        '''Determine total partial charge of molecule'''
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
