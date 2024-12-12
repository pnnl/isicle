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


class XYZGeometryInterface(metaclass=abc.ABCMeta):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'dft')
                and callable(subclass.dft)
                and hasattr(subclass, 'md')
                and callable(subclass.md)
                and hasattr(subclass, 'get_natoms')
                and callable(subclass.get_natoms)
                and hasattr(subclass, 'get_atom_indices')
                and callable(subclass.get_atom_indices)
                and hasattr(subclass, 'get___dict__')
                and callable(subclass.get___dict__)
                and hasattr(subclass, '__copy__')
                and callable(subclass.__copy__)
                and hasattr(subclass, 'to_xyzblock')
                and callable(subclass.to_xyzblock)
                or NotImplemented)

    @abc.abstractmethod
    def dft(self):
        """Optimize geometry using density function theory (DFT) methods."""
        raise NotImplementedError

    @abc.abstractmethod
    def md(self):
        """Optimize geometry using molecule dynamics methods (MD)."""
        raise NotImplementedError

    @abc.abstractmethod
    def get_natoms(self):
        """Count number of atoms"""
        raise NotImplementedError

    @abc.abstractmethod
    def get_atom_indices(self):
        """Extract indices of each atom from the internal geometry."""
        raise NotImplementedError

    @abc.abstractmethod
    def get___dict__(self):
        """Return a copy of this object's attributes dictionary"""
        raise NotImplementedError

    @abc.abstractmethod
    def __copy__(self):
        """Return hard copy of this class instance."""
        raise NotImplementedError

    @abc.abstractmethod
    def to_xyzblock(self):
        """Get XYZ text for this structure."""
        raise NotImplementedError


class GeometryInterface(XYZGeometryInterface):

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'get___dict__')
                and callable(subclass.get___dict__)
                and hasattr(subclass, '__copy__')
                and callable(subclass.__copy__)
                and hasattr(subclass, 'get_total_partial_charge')
                and callable(subclass.get_total_partial_charge)
                and hasattr(subclass, 'dft')
                and callable(subclass.dft)
                and hasattr(subclass, 'md')
                and callable(subclass.md)
                and hasattr(subclass, 'get_natoms')
                and callable(subclass.get_natoms)
                and hasattr(subclass, 'get_atom_indices')
                and callable(subclass.get_atom_indices)
                and hasattr(subclass, 'to_mol')
                and callable(subclass.to_mol)
                and hasattr(subclass, 'to_smiles')
                and callable(subclass.to_smiles)
                and hasattr(subclass, 'to_inchi')
                and callable(subclass.to_inchi)
                and hasattr(subclass, 'to_smarts')
                and callable(subclass.to_smarts)
                and hasattr(subclass, 'to_xyzblock')
                and callable(subclass.to_xyzblock)
                or NotImplemented)

    @abc.abstractmethod
    def to_mol(self, path: str):
        """Returns RDKit Mol object for this Geometry."""
        raise NotImplementedError

    @abc.abstractmethod
    def get_total_partial_charge(self):
        """Determine total partial charge of molecule"""
        raise NotImplementedError

    @abc.abstractmethod
    def to_smiles(self, path: str):
        """Return SMILES representation"""
        raise NotImplementedError

    @abc.abstractmethod
    def to_inchi(self, path: str):
        """Return InChI representation"""
        raise NotImplementedError

    @abc.abstractmethod
    def to_smarts(self, path: str):
        """Return SMARTS representation"""
        raise NotImplementedError


class WrapperInterface(metaclass=abc.ABCMeta):
    """
    Abstract base class for wrapper interfaces.
    """

    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'set_geometry')
                and callable(subclass.set_geometry)
                and hasattr(subclass, 'configure')
                and callable(subclass.configure)
                and hasattr(subclass, 'submit')
                and callable(subclass.submit)
                and hasattr(subclass, 'run')
                and callable(subclass.run)
                and hasattr(subclass, 'finish')
                and callable(subclass.finish)
                or NotImplemented)

    @abc.abstractmethod
    def set_geometry(self):
        """Assign the input geometry."""
        raise NotImplementedError

    @abc.abstractmethod
    def configure(self):
        """Configure the run."""
        raise NotImplementedError

    @abc.abstractmethod
    def submit(self):
        """Configure the run."""
        raise NotImplementedError

    @abc.abstractmethod
    def run(self):
        """Execute/submit the run."""
        raise NotImplementedError

    @abc.abstractmethod
    def finish(self, path: str):
        """Finalize, parse, return result object."""
        raise NotImplementedError
