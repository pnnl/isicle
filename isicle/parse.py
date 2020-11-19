from isicle.interfaces import FileParserInterface
import pandas as pd

class NWChemResult():
    """Organize parsed results from NWChem outputs"""

    def __init__(self):
        self.energy = None  # Dictionary, keys: energy, charges
        self.geometry = None  # String, filename (for now)
        self.shielding = None  # DataFrame
        self.spin = None  # Not set
        self.frequency = None # Dictionary, see function for keys
        self.molden = None  # String, filename (for now)

    def set_energy(self, energy):
        result = {'energy':[energy[0]], 'charges':energy[1]}
        self.energy = energy
        return self.energy

    def set_geometry(self, geometry_filename):
        # TODO: save geometry object instead
        self.geometry = geometry_filename
        return self.geometry

    def set_shielding(self, shielding):

        shield_values, dft_energy, index = shielding

        # TODO: change how this info is stored?
        df = pd.DataFrame(shield_values, columns=['index', 'atom', 'shielding'])
        df['dft_energy'] = energy[-1]
        df['index'] = true_idx
        return self.shielding

    def set_spin(self, spin):
        self.spin = spin
        return self.spin

    def set_frequency(self, frequency):
        '''
        Create dictionary from results and save as attribute.

        Keys: 'natoms', 'lowdinIdx', 'energies', 'enthalpies', 'entropies',
               'capacities', 'preoptTime', 'geomoptTime', 'cpuTime', 'zpe'
        '''

        # Make dictionary with results
        frequency_d = {}
        names = ['natoms', 'lowdinIdx', 'energies', 'enthalpies', 'entropies',
                'capacities', 'preoptTime', 'geomoptTime', 'cpuTime', 'zpe']
        for i, name in enumerate(names):
            frequency_d[name] = frequency[i]

        self.frequency = frequency_d
        return self.frequency

    def set_molden(self, molden_filename):
        # TODO: any processing on file contents?
        self.molden = molden_filename
        return self.molden

    def get_energy(self):
        return self.energy

    def get_geometry(self):
        return self.geometry

    def get_shielding(self):
        return self.shielding

    def get_spin(self):
        return self.spin()

    def get_frequency(self):
        '''
        Return dictionary with frequency-related results

        Keys: 'natoms', 'lowdinIdx', 'energies', 'enthalpies', 'entropies',
               'capacities', 'preoptTime', 'geomoptTime', 'cpuTime', 'zpe'
        '''
        return self.frequency

    def get_molden(self):
        return self.molden


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

        return geoms[-1]

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

        return energy, df.Charge.tolist()

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

        return shield_values, energy[-1], true_idx

    def _parse_spin(self):
        coor_substr='Output coordinates in angstroms'
        cst_substr='Total Shielding Tensor'

        # Extracting Number of Isotopes/Atoms
        # Must sit alone due to pulling number of atoms for later parsing
        natoms = 0
        for i in range(length(self.contents) - 1):  # Minus one since always looking one ahead
           currentline = self.contents[i]
           nextline = self.contents[i + 1]
           if 'property' in currentline and 'SHIELDING' in nextline:
               natoms=int(nextline.split(' ')[-1])

        # Check that natoms was found, exit otherwise
        if natoms == 0:
            print('Number of atoms not found or equal to zero')
            return None

        # Declaring couplings
        coup_freqs = np.zeros((natoms, natoms))

        cst = []
        # Search for spin-spin couplings and populate matrix
        for ii in range(length(self.contents)):
            currentline = self.contents[ii]
            temp = currentline.split(' ')

            # Extracting Isotopes/Atoms & Coordinates
            if coor_substr in currentline:
                # Indexing ii+3+natoms; up to and including
                coor_cellarr = self.contents[ii + 4: ii + 3 + natoms]

            # Extracting Full Chemical Shielding Tensors
            # NOTE: Shielding Tensors only calculated for C&H; must match in
            # handle_parsed.m
            elif cst_substr in currentline:
                cst.append(self.contents[ii + 1])
                cst.append(self.contents[ii + 2])
                cst.append(self.contents[ii + 3])

            # Extracting Coupling Frequencies matrix
            elif temp[0] == 'Atom' and temp[4] == 'Atom':
                col_idx = int(temp[1].replace(':', '')) # Needs to be int
                row_idx = int(temp[5].replace(':', '')) # Needs to be int
                temp_freq = nwctext[ii+41].split(' ')
                if 'Spin-Spin' not in temp_freq[1]:
                    print('Exact line for coupling frequencies not found: Change search parameters in line above.')
                    return None
                temp_freq = float(temp_freq[4])

                # Accounting for symmetry by switching row/column index
                coup_freqs[row_idx, col_idx] += temp_freq
                coup_freqs[col_idx, row_idx] += temp_freq

        # Assert size of coordinate cell array matches number of atoms
        if length(coor_cellarr) != natoms:
            print('Coordinate cell array size does not match number of atoms.')
            return None

        # Ensuring diaganolized zeros
        for ii in range(natoms):
            if coup_freqs[ii, ii]!=0:
                print('Extracted Coupling Frequency incorrect: overwriting to zero.')
                coup_freqs[ii, ii] = 0

        return coup_freqs

    def _parse_frequency(self):

        # Init
        indices = []
        lowdinIdx = []
        energies = []
        zpe = []
        enthalpies = []
        entropies = []
        capacities = []
        preoptTime = 0
        geomoptTime = 0
        freqTime = 0
        cpuTime = 0
        wallTime = 0
        ready = False
        opt = False
        freq = False

        for i, line in enumerate(self.contents):

            # ?
            if 'No.' in line and len(indices) == 0:
                indices.append(i + 2)  # 0
            elif 'Atomic Mass' in line and len(indices) == 1:
                indices.append(i - 1)  # 1
                indices.append(i + 3)  # 2
            elif 'Effective nuclear repulsion energy' in line and len(indices) == 3:
                indices.append(i - 2)  # 3
            elif 'Optimization converged' in line:  # lowdin population
                ready = True
            elif 'Failed to converge in maximum number of steps or available time' in line:
                ready=True
            elif 'Output coordinates in angstroms' in line:  # Shell charges
                lowdinIdx.append(i + 4)
                ready = False

            # Include? Commented or from past files
            #elif ready is True:
                #lowdinIdx.append(i + 2)
                #ready = False
            elif 'Shell Charges' in line and ready is True: ##Shell Charges
                lowdinIdx.append(i + 2)
                ready = False
            elif 'Lowdin Population Analysis' in line:
                ready = True

            # Get values
            if 'Total DFT energy' in line:
                energies.append(float(line.rstrip().split('=')[-1]))

            if 'Zero-Point correction to Energy' in line:
                zpe.append(line.rstrip().split('=')[-1])

            if 'Thermal correction to Enthalpy' in line:
                enthalpies.append(line.rstrip().split('=')[-1])

            if 'Total Entropy' in line:
                entropies.append(line.rstrip().split('=')[-1])

            if 'constant volume heat capacity' in line:
                capacities.append(line.rstrip().split('=')[-1])

            # Check for optimization and frequency calcs
            if 'NWChem Geometry Optimization' in line:
                opt = True
            elif 'NWChem Nuclear Hessian and Frequency Analysis' in line:
                freq = True

            # Get timing
            if 'Total iterative time' in line and opt is False:
                preoptTime += float(line.rstrip().split('=')[1].split('s')[0])
            elif 'Total iterative time' in line and opt is True and freq is False:
                geomoptTime += float(line.rstrip().split('=')[1].split('s')[0])
            elif 'Total iterative time' in line and freq is True:
                freqTime += float(line.rstrip().split('=')[1].split('s')[0])

            if 'Total times' in line:
                cpuTime = float(line.rstrip().split(':')[1].split('s')[0])
                wallTime = float(line.rstrip().split(':')[2].split('s')[0])
                freqTime = (cpuTime - geomoptTime - preoptTime)

        natoms = int(self.contents[indices[1] - 1].split()[0])

        return natoms, lowdinIdx, energies, enthalpies, entropies, capacities, \
               preoptTime, geomoptTime, cpuTime, zpe

    def _parse_mulliken(self, path):

        # search = splitext(path)[0]
        # m = glob.glob(search + '*.molden')
        #
        # if len(m) != 1:
        #     raise IOError('Incorrect number of molden files found.')
        #
        # return m

        return None

    # TODO: what should default to_parse be?
    def parse(self, to_parse=['geometry', 'energy', 'shielding', 'spin', 'frequency'],
              geom_path=self.path):
        """Extract relevant information from data"""

        result = NWChemResult()

        if 'geometry' in to_parse:

            try:
                geometry_filename = _parse_geometry_filename(geom_path)
                result.set_geometry(geometry_filename)  # Store as filename

            except IndexError:
                pass

        if 'energy' in to_parse:

            try:
                energy = _parse_energy()
                result.set_energy(energy)  # Stored as dictionary
            except IndexError:
                pass

        if 'shielding' in to_parse:
            try:
                shielding = _parse_shielding()
                result.set_shielding(shielding)  # Stored as dictionary
            except IndexError:
                pass

        if 'spin' in to_parse:  # N2S
            try:
                spin = _parse_spin()
                result.set_spin(spin)
            except IndexError:
                pass

        if 'frequency' in to_parse:
            try:
                frequency = _parse_frequency()
                result.set_frequency(frequency)
            except IndexError:
                pass

        if 'mulliken' in to_parse:
            try:
                mulliken = _parse_mulliken()
                result.set_mulliken(mulliken)
            except IndexError:
                pass

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
