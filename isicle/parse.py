from isicle.interfaces import FileParserInterface
import pandas as pd


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
    def __init__(self):
        self.contents = None
        self.result = None

    def load(self, path: str):
        """Load in the data file"""
        with open(path, 'r') as f:
            self.contents = f.readlines()

        return self.contents

    def parse(self):
        """Extract relevant information from data"""
        done = False
        for line in self.contents:
            # if "average (second order) TM mobility" in line:
            #     m_mn = float(line.split('=')[-1])
            if "average TM cross section" in line:
                ccs_mn = float(line.split('=')[-1])
            elif "standard deviation TM cross section" in line:
                ccs_std = float(line.split('=')[-1])
            elif 'standard deviation (percent)' in line:
                done = True
        if done is True:
            self.result = {'ccs': [ccs_mn], 'std': [ccs_std]}

        return self.result

    def save(self, path: str, sep='\t'):
        """Write parsed object to file"""
        pd.DataFrame(self.result).to_csv(path, sep=sep, index=False)


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
