import abc


class FileParserInterface(metaclass=abc.ABCMeta):
    """
    Abstract base class for file parser interface. All file parsers
    conform to this definition.
    """

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
        """Load in the data file"""
        raise NotImplementedError

    @abc.abstractmethod
    def parse(self):
        """Extract relevant information from data"""
        raise NotImplementedError

    @abc.abstractmethod
    def save(self, path: str):
        """Write parsed object to file"""
        raise NotImplementedError


class MolecularStringInterface(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'load')
                and callable(subclass.load)
                and hasattr(subclass, 'desalt')
                and callable(subclass.desalt)
                and hasattr(subclass, 'neutralize')
                and callable(subclass.neutralize)
                and hasattr(subclass, 'tautomerize')
                and callable(subclass.tautomerize)
                or NotImplemented)

    @abc.abstractmethod
    def load(self, path: str):
        """Load in the molecular data file"""
        raise NotImplementedError

    @abc.abstractmethod
    def desalt(self):
        """Desalt molecular object"""
        raise NotImplementedError

    @abc.abstractmethod
    def neutralize(self):
        """Neutralize molecular object"""
        raise NotImplementedError

    @abc.abstractmethod
    def tautomerize(self):
        """Tautomerize molecular object"""
        raise NotImplementedError


class GeometryInterface(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'load')
                and callable(subclass.load)
                and hasattr(subclass, 'convert2D')
                and callable(subclass.convert2D)
                and hasattr(subclass, 'convert3D')
                and callable(subclass.convert3D)
                and hasattr(subclass, 'save')
                and callable(subclass.save)
                or NotImplemented)

    @abc.abstractmethod
    def load(self, path: str):
        """Load in SMILES string"""
        raise NotImplementedError

    @abc.abstractmethod
    def convert2D(self):
        """Convert into 2D geometry using RDKit"""
        raise NotImplementedError

    @abc.abstractmethod
    def convert3D(self):
        """Convert 2D geometry into 3D geometry using RDKit"""
        raise NotImplementedError

    @abc.abstractmethod
    def save(self, path: str):
        """Write 3D molecule to file (xyz and mol/mol2)"""
        raise NotImplementedError

    @abc.abstractmethod
    def natoms(self):
        """Count number of atoms"""
        raise NotImplementedError

    @abc.abstractmethod
    def total_partial_charge(self):
        raise NotImplementedError

    @abc.abstractmethod
    def pop_atom(path, output, atom='Na'):
        """Removes Na atoms"""
        raise NotImplementedError

    @abc.abstractmethod
    def push_atom(path, output, idx, content):
        raise NotImplementedError

class AdductInterface(metaclass=abc.ABCMeta):
    
