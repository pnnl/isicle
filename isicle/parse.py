#from isicle.interfaces import FileParserInterface
import pandas as pd

class NWChemParser(FileParserInterface):
    """Extract text from an NWChem simulation output file."""

    def load(self, path: str):
        """Load in the data file"""
        raise NotImplementedError

    def parse(self, to_parse=['geometry', 'energy', 'shielding', 'spin']):
        """Extract relevant information from data"""
        done = True  # TODO: check for error file (need example)
        for line in self.contents:
            pass
        raise NotImplementedError

    def save(self, path: str):
        """Write parsed object to file"""
        raise NotImplementedError


class ImpactParser(FileParserInterface):
    """Extract text from an Impact mobility calculation output file."""

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

        # Check CCS results == 1
        count = 0
        for line in self.contents:
            l = line.split(' ')
            if 'CCS' in l[0]:
                count += 1
        if count != 1:
            return self.result

        # Assume values in second line
        l = self.contents[1].split(' ')
        l = [x for x in l if len(x) > 0]

        # Pull values of interest - may be error prone
        values = []
        try:
            values.append(float(l[-5]))
            values.append(float(l[-3][:-1]))
            values.append(float(l[-2]))
            values.append(int(l[-1]))
        except (ValueError, IndexError) as e:
            print('Could not parse file: ', e)
            return None

        # Add to dictionary to return
        result = {}
        keys = ['CCS_PA', 'SEM_rel', 'CCS_TJM', 'n_iter']
        for key, val in zip(keys, values):
            result[key] = [val]

        # Save and return results
        self.result = result
        return result # TODO: return CCS?

    def save(self, path: str, sep='\t'):
        """Write parsed object to file"""
        pd.DataFrame(self.result).to_csv(path, sep=sep, index=False)
        return


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
        return


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
