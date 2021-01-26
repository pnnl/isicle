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
        return (hasattr(subclass, '_to_2D')
                and callable(subclass._to_2D)
                and hasattr(subclass, 'to_2D')
                and callable(subclass.to_2D)
                and hasattr(subclass, '_to_3D')
                and callable(subclass._to_3D)
                and hasattr(subclass, 'to_3D')
                and callable(subclass.to_3D)
                and hasattr(subclass, 'natoms')
                and callable(subclass.natoms)
                and hasattr(subclass, 'total_partial_charge')
                and callable(subclass.total_partial_charge)
                and hasattr(subclass, 'save_pickle')
                and callable(subclass.save_pickle)
                and hasattr(subclass, 'save_mol')
                and callable(subclass.save_mol)
                and hasattr(subclass, 'copy')
                and callable(subclass.copy)
                or NotImplemented)

    @abc.abstractmethod
    def _to_2D(self, path: str):
        """Load in SMILES string"""
        raise NotImplementedError

    @abc.abstractmethod
    def to_2D(self):
        """Convert into 2D geometry using RDKit"""
        raise NotImplementedError

    @abc.abstractmethod
    def _to_3D(self):
        """Convert 2D geometry into 3D geometry using RDKit"""
        raise NotImplementedError

    @abc.abstractmethod
    def to_3D(self):
        """Convert 2D geometry into 3D geometry using RDKit"""
        raise NotImplementedError

    @abc.abstractmethod
    def save_pickle(self, path: str):
        raise NotImplementedError

    @abc.abstractmethod
    def save_mol(self, path: str):
        """Write 3D molecule to file (xyz and mol/mol2)"""
        raise NotImplementedError

    @abc.abstractmethod
    def natoms(self):
        """Count number of atoms"""
        raise NotImplementedError

    @abc.abstractmethod
    def total_partial_charge(self):
        raise NotImplementedError

class AdductInterface(metaclass=abc.ABCMeta):
    
