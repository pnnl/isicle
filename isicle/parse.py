#from isicle.interfaces import FileParserInterface
import pandas as pd

class NWChemParser(FileParserInterface):
    """Extract text from an NWChem simulation output file."""

    def __init__(self):
        self.contents = None
        self.result = None
        self.path = None

    def load(self, path: str):
        """Load in the data file"""
        with open(path, 'r') as f:
            self.contents = f.readlines()
        self.path = path
        return self.contents

    # TODO: use Amy's load_geom to generate object
    def _parse_geometry(self, path):
        search = splitext(path)[0]
        geoms = glob.glob(search + '*.xyz')

        if len(geoms) < 1:
            raise IOError('No geometry files found.')

        geoms.sort()

        # Create final dict to return
        result = {'geometry':geoms[-1]]}

        return result

    def _parse_energy(self):

        energy = []
        charges = []
        ready = False
        for line in self.contents:
            if 'Total DFT energy' in line:
                energy.append(float(line.split()[-1]))

            elif 'Atom       Charge   Shell Charges' in line:
                ready = True
                charges = []
            elif ready is True and line.strip() in ['', 'Line search:']:
                ready = False
            elif ready is True:
                charges.append(line)

        # grab last energy
        energy = energy[-1]

        # process charge information
        df = pd.DataFrame([x.split()[0:4] for x in charges[1:]],
                          columns=['idx', 'Atom', 'Number', 'Charge'])

        df.Number = df.Number.astype('int')
        df.Charge = df.Number - df.Charge.astype('float')

        # Create final dict to return
        result = {'energy':[energy], 'charges':df.Charge.tolist()}

        return result

    def _parse_shielding(self):

        # Init
        energy = []
        shield_values = []
        ready = False

        # TODO: maybe this?
        shield_idxs = []
        shield_atoms = []
        shields = []

        for line in self.contents:
            if "Total DFT energy" in line:  # TODO: should this go elsewhere?
                energy.append(float(line.split()[-1]))
            elif "Atom:" in line:
                idx = int(line.split()[1])
                atom = line.split()[2]
                ready = True
            elif "isotropic" in line and ready is True:
                shield = float(line.split()[-1])
                shield_values.append([idx, atom, shield])
                ready = False
            elif 'SHIELDING' in line:
                true_idx = [int(x) for x in line.split()[2:]]

        #df = pd.DataFrame(shield_values, columns=['index', 'atom', 'shielding'])
        #df['dft_energy'] = energy[-1]
        #df['index'] = true_idx

        # Creare final dict to return
        result = {}
        # TODO: what is needed?
        result['shields'] = shields
        result['shield_atoms'] = shield_atoms
        result['shield_idxs'] = shield_idxs
        result['shielding_dft_energy'] = energy[-1]
        result['shielding_index'] = true_idx

        return result

    def _parse_spin(self):
        return

    def parse(self, to_parse=['geometry', 'energy', 'shielding', 'spin'], geom_path=self.path):
        """Extract relevant information from data"""

        result = {}

        if 'geometry' in to_parse:
            result.update(_parse_geometry(geom_path))

        if 'energy' in to_parse:
            result.update(_parse_energy())

        if 'shielding' in to_parse:
            result.update(_parse_shielding())

        if 'spin' in to_parse:  # N2S
            result.update(_parse_spin())

        # Add frequency parser

        return result

    def save(self, path: str):
        """Write parsed object to file"""
        pd.DataFrame(self.result).to_csv(path, sep=sep, index=False)
        return


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
