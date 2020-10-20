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
