from isicle.interfaces import FileParserInterface
import pandas as pd
from os.path import splitext
import glob
import pickle
import numpy as np


class NWChemResult():
    '''Organize parsed results from NWChem outputs'''

    def __init__(self):
        self.energy = None  # Dictionary, keys: energy, charges
        self.geometry = None  # String, filename (for now)
        self.shielding = None  # Dictionary
        self.spin = None  # Not set
        self.frequency = None  # Dictionary, see function for keys
        self.molden = None  # String, filename (for now)
        self.meta = None  # Dictionary, see function for keys

    def set_energy(self, energy):
        result = {'energy': [energy[0]]}
        self.energy = result
        return self.energy

    def set_geometry(self, geometry_filename):
        # TODO: save xyz block as well
        self.geometry = geometry_filename
        return self.geometry

    def set_shielding(self, shielding):
        result = {'index': shielding[0], 'atom': shielding[1],
                  'shielding': shielding[2]}
        self.shielding = result
        return self.shielding

    def set_spin(self, spin):
        # TODO
        self.spin = spin
        return self.spin

    def set_frequency(self, frequency):
        self.frequency = frequency
        return self.frequency

    def set_timing(self, timing):
        result = {'single point': timing[0], 'geometry optimization': timing[1],
                  'frequency': timing[2], 'total': timing[3]}
        self.timing = result
        return self.timing

    def set_charge(self,charge):
        result = {'charge': charge}
        self.charge = result
        return self.charge

    def set_meta(self, meta):
        '''
        Create dictionary from results and save as attribute.

        Keys: 'natoms', 'lowdinIdx', 'energies', 'enthalpies', 'entropies',
               'capacities', 'preoptTime', 'geomoptTime', 'cpuTime', 'zpe'
        '''

        # Make dictionary with results
        meta_d = {}
        names = ['natoms', 'lowdinIdx', 'energies', 'enthalpies', 'entropies',
                 'capacities', 'preoptTime', 'geomoptTime', 'cpuTime', 'zpe']
        for i, name in enumerate(names):
            meta_d[name] = meta[i]

        self.meta = meta_d
        return self.meta

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
        return self.spin

    def get_frequency(self):
        return self.frequency

    def get_timing(self):
        return self.timing

    def get_charge(self):
        return self.charge

    def get_molden(self):
        return self.molden

    def save(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        return

    def load(self, path):
        '''
        Load saved class data to this (overwrites all variables)
        '''

        # Load existing file
        with open(path, 'rb') as f:
            saved_result = pickle.load(f)

        # Overwrite the variables in this object
        self.geometry = saved_result.get_geometry()
        self.energy = saved_result.get_energy()
        self.shielding = saved_result.get_shielding()
        self.spin = saved_result.get_spin()
        self.frequency = saved_result.get_frequency()
        self.molden = saved_result.get_molden()
        self.timing = saved_result.get_timing()
        self.charge = saved_result.get_charge()
        return

    def to_dict(self):
        d = {}

        d['geometry'] = self.geometry
        d['energy'] = self.energy
        d['shielding'] = self.shielding
        d['spin'] = self.spin
        d['frequency'] = self.frequency
        d['molden'] = self.molden
        d['timing'] = self.timing
        d['charge'] = self.charge

        return d


class NWChemParser(FileParserInterface):
    '''Extract text from an NWChem simulation output file.'''

    def __init__(self):
        self.contents = None
        self.result = None
        self.path = None

    def load(self, path: str):
        '''Load in the data file'''
        with open(path, 'r') as f:
            self.contents = f.readlines()
        self.path = path
        return self.contents

    def _parse_geometry_filename(self, path):
        '''Grab path to .xyz file or generate .xyz file from *.out file '''
        search = splitext(path)[0]
        geoms = glob.glob(search + '*.xyz')
        coor_substr = 'Output coordinates in angstroms'

        # Extracting Atoms & Coordinates
        ii = [i for i in range(len(self.contents)) if coor_substr in self.contents[i]]
        ii.sort()

        coord = ''
        g = ii[-1]+4
        natoms = 0
        while g <= len(self.contents)-1:
            if self.contents[g] != ' \n':
                line = self.contents[g].split()
                xyz_line = line[1] + '\t' + line[3] + '\t' + line[4] + '\t' +line[5] +'\n'
                coord += xyz_line
                natoms += 1

            else:
                break
            g+=1

        coord = str(natoms) + '\n\n' + coord 
        name = search + '.xyz'
        xyz_file = open(name, 'w')
        f = xyz_file.write(coord)
        xyz_file.close()
        
        return name

    def _parse_energy(self):

        #TO DO: Add Initial energy and final energy if different

        # Init
        energy = None

        # Cycle through file
        for line in self.contents:
            if 'Total DFT energy' in line:
                # Overwrite last saved energy
                energy = line.split()[-1]

        return energy, None

    def _parse_shielding(self):

        # Init
        ready = False
        shield_idxs = []
        shield_atoms = []
        shields = []

        for line in self.contents:
            if "Atom:" in line:
                idx = int(line.split()[1])
                atom = line.split()[2]
                ready = True
            elif "isotropic" in line and ready is True:
                shield = float(line.split()[-1])
                shield_idxs.append(idx)
                shield_atoms.append(atom)
                shields.append(shield)
                ready = False

        return shield_idxs, shield_atoms, shields

    def _parse_spin(self):
        # TO DO: Add g-factors

        coor_substr = 'Output coordinates in angstroms'
        cst_substr = 'Total Shielding Tensor'

        # Extracting Number of Isotopes/Atoms
        # Must sit alone due to pulling number of atoms for later parsing
        natoms = 0
        # Minus one since always looking one ahead
        for i in range(len(self.contents) - 1):
            currentline = self.contents[i]
            nextline = self.contents[i + 1]
            if 'property' in currentline and 'SHIELDING' in nextline:
                natoms = int(nextline.split(' ')[-1])

        # Check that natoms was found, exit otherwise
        if natoms == 0:
            print('Number of atoms not found or equal to zero')
            return None

        # Declaring couplings
        coup_freqs = np.zeros((natoms, natoms))
        coup_pairs = []
        coup = []
        ready = False

        for line in self.contents:
            if "Atom  " in line:
                line = line.split()
                idx1 = int((line[1].split(":"))[0]) - 1
                idx2 = int((line[5].split(":"))[0]) - 1
                ready=True
            elif "Isotropic Spin-Spin Coupling =" in line and ready is True:
                coup = float(line.split()[4])
                coup_freqs[idx1][idx2] += coup
                coup_freqs[idx2][idx1] += coup

        # Ensuring diaganolized zeros
        for ii in range(natoms):
            if coup_freqs[ii, ii] != 0:
                print('Extracted Coupling Frequency incorrect: overwriting to zero.')
                coup_freqs[ii, ii] = 0

        return coup_freqs

    def _parse_frequency(self):
        # TO DO: Add Thermo information (zpe, enthalpy, entropy, capacity, rotational constants) 
        #        into the energy global property
        # TO DO: Add freq intensities
        # TO DO: Add rotational/translational/vibrational Cv and entropy
        energies = []
        zpe = []
        enthalpies = []
        entropies = []
        capacities = []
        temp = []
        scaling = []
        natoms = None

        for i, line in enumerate(self.contents):
            if ('Geometry' in line) and (natoms is None):
                atom_start = i + 7
            if ('Atomic Mass' in line) and (natoms is None):
                atom_stop = i - 2
                natoms = atom_stop - atom_start + 1
            if 'Normal Eigenvalue' in line:
                freq_start = i + 3
                freq_stop = i + 2 + 3 * natoms

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

        return np.array([float(x.split()[1])
                         for x in self.contents[freq_start:freq_stop + 1]]), \
               energies, enthalpies, entropies, capacities, zpe

    def _parse_charge(self):
        # TO DO: Parse molecular charge and atomic charges
        # TO DO: Add type of charge
        # TO DO: Multiple instances of charge analysis seen (two Mulliken and one Lowdin, difference?)
        charges = []
        ready = False

        for line in self.contents:

            # Load charges from table
            if 'Atom       Charge   Shell Charges' in line:
                # Table header found. Overwrite anything saved previously
                ready = True
                charges = []
            elif ready is True and line.strip() in ['', 'Line search:']:
                # Table end found
                ready = False
            elif ready is True:
                # Still reading from charges table
                charges.append(line)

            # Include? Commented or from past files
            # elif ready is True:
            #     lowdinIdx.append(i + 2)
            #     ready = False
            elif 'Shell Charges' in line and ready is True:  # Shell Charges
                lowdinIdx.append(i + 2)
                ready = False
            elif 'Lowdin Population Analysis' in line:
                ready = True

        # Process table if one was found
        if len(charges) > 0:

            # Remove blank line in charges (table edge)
            charges = charges[1:]

            # Process charge information
            df = pd.DataFrame([x.split()[0:4] for x in charges],
                              columns=['idx', 'Atom', 'Number', 'Charge'])
            df.Number = df.Number.astype('int')
            df.Charge = df.Number - df.Charge.astype('float')

            return energy, df.Charge.tolist()

        return None

    def _parse_timing(self):

        # Init
        indices = []
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

        return preoptTime, geomoptTime, freqTime, cpuTime

    def _parse_molden(self, path):

        search = splitext(path)[0]
        m = glob.glob(search + '*.molden')

        if len(m) != 1:
            return None

        return m[0]

    def _parse_protocol(self):
        raise NotImplementedError

    # TODO: what should default to_parse be?
    def parse(self, to_parse=['geometry', 'energy'],
              geom_path=None, molden_path=None):
        '''Extract relevant information from data'''

        # Check that the file is valid first
        if len(self.contents) == 0:
            raise RuntimeError('No contents to parse: {}'.format(self.path))
        if 'Total times  cpu' not in self.contents[-1]:
            raise RuntimeError('Incomplete NWChem run: {}'.format(self.path))

        # Initialize result object to store info
        result = NWChemResult()

        if 'geometry' in to_parse:

            try:
                if geom_path is None:
                    geom_path = self.path
                geometry_filename = self._parse_geometry_filename(geom_path)
                result.set_geometry(geometry_filename)  # Store as filename

            except IndexError:
                pass

        if 'energy' in to_parse:

            try:
                energy = self._parse_energy()
                result.set_energy(energy)  # Stored as dictionary
            except IndexError:
                pass

        if 'shielding' in to_parse:
            try:
                shielding = self._parse_shielding()
                result.set_shielding(shielding)  # Stored as dictionary
            except UnboundLocalError:  # Must be no shielding info
                pass

        if 'spin' in to_parse:  # N2S
            try:
                spin = self._parse_spin()
                result.set_spin(spin)
            except IndexError:
                pass

        if 'frequency' in to_parse:
            try:
                frequency = self._parse_frequency()
                result.set_frequency(frequency)
            except IndexError:
                pass

        if 'molden' in to_parse:
            try:
                if molden_path is None:
                    molden_path = self.path
                molden_filename = self._parse_molden(molden_path)
                result.set_molden(molden_filename)
            except IndexError:
                pass

        if 'charge' in to_parse:
            try:
                charge = self._parse_charge()
                result.set_charge(charge)
            except IndexError:
                pass

        if 'timing' in to_parse:
            try:
                timing = self._parse_timing()
                result.set_timing(timing)
            except IndexError:
                pass

        self.result = result
        return result

    def save(self, path: str):
        '''Write parsed object to file'''
        self.result.save(path)
        return


class ImpactParser(FileParserInterface):
    '''Extract text from an Impact mobility calculation output file.'''

    def __init__(self):
        self.contents = None
        self.result = None

    def load(self, path: str):
        '''Load in the data file'''
        with open(path, 'r') as f:
            self.contents = f.readlines()

        return self.contents

    def parse(self):
        '''Extract relevant information from data'''

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
        return result  # TODO: return CCS?

    def save(self, path: str, sep='\t'):
        '''Write parsed object to file'''
        pd.DataFrame(self.result).to_csv(path, sep=sep, index=False)
        return


class MobcalParser(FileParserInterface):
    '''Extract text from a MOBCAL mobility calculation output file.'''

    def __init__(self):
        self.contents = None
        self.result = None

    def load(self, path: str):
        '''Load in the data file'''
        with open(path, 'r') as f:
            self.contents = f.readlines()

        return self.contents

    def parse(self):
        '''Extract relevant information from data'''
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
        '''Write parsed object to file'''
        pd.DataFrame(self.result).to_csv(path, sep=sep, index=False)
        return


class SanderParser(FileParserInterface):
    '''Extract text from an Sander simulated annealing simulation output file.'''

    def load(self, path: str):
        '''Load in the data file'''
        raise NotImplementedError

    def parse(self):
        '''Extract relevant information from data'''
        raise NotImplementedError

    def save(self, path: str):
        '''Write parsed object to file'''
        raise NotImplementedError
