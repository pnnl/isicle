from isicle.interfaces import FileParserInterface


class NWChemParser(FileParserInterface):
    """Extract text from an NWChem simulation output file."""

    def load(self, path: str):
        """Load in the data file"""
        raise NotImplementedError

    def parse(self, to_parse=['geometry', 'energy', 'shielding', 'spin']):
        """Extract relevant information from data"""
        raise NotImplementedError

    def save(self, path: str):
        """Write parsed object to file"""
        raise NotImplementedError


class ImpactParser(FileParserInterface):
    """Extract text from an Impact mobility calculation output file."""

    def load(self, path: str):
        """Load in the data file"""
        raise NotImplementedError

    def parse(self):
        """Extract relevant information from data"""
        raise NotImplementedError

    def save(self, path: str):
        """Write parsed object to file"""
        raise NotImplementedError


class MobcalParser(FileParserInterface):
    """Extract text from a MOBCAL mobility calculation output file."""

    def load(self, path: str):
        """Load in the data file"""
        raise NotImplementedError

    def parse(self):
        """Extract relevant information from data"""
        raise NotImplementedError

    def save(self, path: str):
        """Write parsed object to file"""
        raise NotImplementedError


class SanderParser(FileParserInterface):
    """Extract text from an Sander simulated annealing simulation output file."""

    def load(self, path: str):
        """Load in the data file"""
        raise NotImplementedError

    def parse(self):
        """Extract relevant information from data"""
        raise NotImplementedError

    def save(self, path: str):
        """Write parsed object to file"""
        raise NotImplementedError
