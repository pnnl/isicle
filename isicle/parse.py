from isicle.interfaces import FileParserInterface
import pandas as pd
from os.path import splitext
import glob
import os
import pickle
import numpy as np
from openbabel import pybel
import isicle


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

    def _parse_geometry(self):
        search = os.path.dirname(self.path)
        geoms = sorted(glob.glob(os.path.join(search, '*.xyz')))
        
        if len(geoms) > 0:
            return isicle.geometry.load(geoms[-1])

        raise Exception

    def _parse_energy(self):

        # TO DO: Add Initial energy and final energy if different

        # Init
        energy = None

        # Cycle through file
        for line in self.contents:
            if 'Total DFT energy' in line:
                # Overwrite last saved energy
                energy = float(line.split()[-1])

        return energy

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

        if len(shields) > 1:
            return {'index': shield_idxs,
                    'atom': shield_atoms,
                    'shielding': shields}
        
        raise Exception

    def _parse_spin(self):
        # TO DO: Add g-factors

        # Declaring couplings
        coup_pairs = []
        coup = []
        index = []
        g_factor = []
        ready = False

        for line in self.contents:
            if "Atom  " in line:
                line = line.split()
                idx1 = int((line[1].split(":"))[0])
                idx2 = int((line[5].split(":"))[0])
                ready = True
            elif "Isotropic Spin-Spin Coupling =" in line and ready is True:
                coupling = float(line.split()[4])
                coup_pairs.append([idx1, idx2])
                coup.append(coup)
                ready = False
            elif "Respective Nuclear g-factors:" in line:
                line = line.split()
                if idx1 not in index:
                    index.append(idx1)
                    g = float(line[3])
                    g_factor.append(g)
                if idx2 not in index:
                    index.append(idx2)
                    g = float(line[5])
                    g_factor.append(g)

        return {'pair indices': coup_pairs, 'spin couplings': coup,
                'index': index, 'g-tensors': g_factor}

    def _parse_frequency(self):
        # TO DO: Add freq intensities
        # TO DO: Add rotational/translational/vibrational Cv and entropy
        freq = None
        zpe = None
        enthalpies = None
        entropies = None
        capacities = None
        temp = None
        scaling = None
        natoms = None
        has_frequency = None

        for i, line in enumerate(self.contents):
            if ('geometry' in line) and (natoms is None):
                atom_start = i + 7
            if ('Atomic Mass' in line) and (natoms is None):
                atom_stop = i - 2
                natoms = atom_stop - atom_start + 1
            if 'Normal Eigenvalue' in line:
                has_frequency = True
                freq_start = i + 3
                freq_stop = i + 2 + 3 * natoms

            # Get values
            if 'Zero-Point correction to Energy' in line:
                zpe = line.rstrip().split('=')[-1]

            if 'Thermal correction to Enthalpy' in line:
                enthalpies = line.rstrip().split('=')[-1]

            if 'Total Entropy' in line:
                entropies = line.rstrip().split('=')[-1]

            if 'constant volume heat capacity' in line:
                capacities = line.rstrip().split('=    ')[-1]

        if has_frequency is True:
            freq = np.array([float(x.split()[1]) for x in self.contents[freq_start:freq_stop + 1]])
            
            return {'frequencies': freq, 'correction to enthalpy': enthalpies,
                    'total entropy': entropies, 'constant volume heat capacity': capacities,
                    'zero-point correction': zpe}
        
        raise Exception

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
            # elif 'Shell Charges' in line and ready is True:  # Shell Charges
            #     lowdinIdx.append(i + 2)
            #     ready = False
            # elif 'Lowdin Population Analysis' in line:
            #     ready = True

        # Process table if one was found
        if len(charges) > 0:
            # return charges

            # Remove blank line in charges (table edge)
            charges = charges[1:]

            # Process charge information
            df = pd.DataFrame([x.split()[0:4] for x in charges],
                              columns=['idx', 'Atom', 'Number', 'Charge'])
            df.Number = df.Number.astype('int')
            df.Charge = df.Number - df.Charge.astype('float')

            # TODO: is this how we want to return?
            return df.Charge.values

        raise Exception

    def _parse_timing(self):

        # Init
        indices = []
        preoptTime = 0
        geomoptTime = 0
        freqTime = 0
        cpuTime = 0
        # wallTime = 0
        # ready = False
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
            if 'NWChem geometry Optimization' in line:
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
                # wallTime = float(line.rstrip().split(':')[2].split('s')[0])
                freqTime = (cpuTime - geomoptTime - preoptTime)

        # natoms = int(self.contents[indices[1] - 1].split()[0])

        return {'single point': preoptTime, 'geometry optimization': geomoptTime,
                'frequency': freqTime, 'total': cpuTime}

    def _parse_molden(self):

        search = splitext(self.path)[0]
        m = glob.glob(search + '*.molden')

        if len(m) != 1:
            raise Exception

        return m[0]

    def _parse_protocol(self):

        '''Parse out dft protocol'''
        functional = []
        basis_set = []
        solvation = []
        tasks = []
        basis = None
        func = None
        solvent = None

        for line in self.contents:

            if '* library' in line:
                basis = line.split()[-1]
            if ' xc ' in line:
                func = line.split(' xc ')[-1].strip()
            if 'solvent ' in line:
                solvent = line.split()[-1]
            if 'task dft optimize' in line:
                tasks.append('optimize')
                basis_set.append(basis)
                functional.append(func)
                solvation.append(solvent)
            if 'SHIELDING' in line:
                tasks.append('shielding')
                basis_set.append(basis)
                functional.append(func)
                solvation.append(solvent)
            if 'SPINSPIN' in line:
                tasks.append('spin')
                basis_set.append(basis)
                functional.append(func)
                solvation.append(solvent)
            if 'freq ' in line:
                tasks.append('frequency')
                basis_set.append(basis)
                functional.append(func)
                solvation.append(solvent)

        return {'functional': functional, 'basis set': basis_set,
                  'solvation': solvation, 'tasks': tasks}

    def parse(self, to_parse=['geometry', 'energy']):
        '''
        Extract relevant information from NWChem output

        Parameters
        ----------
        to_parse : list of str
            geometry, energy, shielding, spin, frequency, molden, charge, timing 
        '''

        # Check that the file is valid first
        if len(self.contents) == 0:
            raise RuntimeError('No contents to parse: {}'.format(self.path))
        if 'Total times  cpu' not in self.contents[-1]:
            raise RuntimeError('Incomplete NWChem run: {}'.format(self.path))

        # Initialize result object to store info
        result = {}

        try:
            result['protocol'] = self._parse_protocol()
        except:
            pass

        if 'geometry' in to_parse:
            try:
                result['geom'] = self._parse_geometry()

            except:
                pass

        if 'energy' in to_parse:
            try:
                result['energy'] = self._parse_energy()
            except:
                pass

        if 'shielding' in to_parse:
            try:
                result['shielding'] = self._parse_shielding()
            except:  # Must be no shielding info
                pass

        if 'spin' in to_parse:  # N2S
            try:
                result['spin'] = self._parse_spin()
            except:
                pass

        if 'frequency' in to_parse:
            try:
                result['frequency'] = self._parse_frequency()
            except:
                pass

        if 'molden' in to_parse:
            try:
                 result['molden'] = self._parse_molden()
            except:
                pass

        if 'charge' in to_parse:
            try:
                result['charge'] = self._parse_charge()
            except:
                pass

        if 'timing' in to_parse:
            try:
                result['timing'] = self._parse_timing()
            except:
                pass

        return result

    def save(self, path):
        with open(path, 'wb') as f:
            pickle.dump(self, f)
        return


class ImpactParser(FileParserInterface):
    '''Extract text from an Impact mobility calculation output file.'''

    def __init__(self):
        self.contents = None
        self.result = None

    def load(self, path: str):
        '''Load in the data file'''
        with open(path, 'rb') as f:
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
        with open(path, 'rb') as f:
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


class XTBResult():

    def __init__(self):
        self.energy = None  # Dictionary, keys: energy, charges
        self.geom = None  # geometry object or array of geometry obejects
        self.timing = None  # Dictionary, see function for keys
        self.protocol = None  # Dictionary

    def set_energy(self, energy):
        result = energy
        self.energy = result
        return self.energy

    def set_geometry(self, geom):
        self.geom = geom
        return self.geom

    def set_timing(self, timing):
        result = timing
        self.timing = result
        return self.timing

    def set_protocol(self, protocol):
        result = protocol
        self.protocol = result
        return self.protocol

    def get_energy(self):
        return self.energy

    def get_geometry(self):
        return self.geom

    def get_timing(self):
        return self.timing

    def get_protocol(self):
        return self.protocol

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
        self.geom = saved_result.get_geometry()
        self.energy = saved_result.get_energy()
        self.timing = saved_result.get_timing()
        self.protocol = saved_result.get_protocol()
        return

    def to_dict(self):
        d = {}

        d['geometry'] = self.geom
        d['energy'] = self.energy
        d['timing'] = self.timing
        d['protocol'] = self.protocol

        return d


class XTBParser(FileParserInterface):
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

    def _crest_energy(self):

        relative_energy = []
        total_energy = []
        population = []

        ready = False
        h = 0
        while h <= len(self.contents)-1:
            if 'Erel/kcal     Etot      weight/tot conformer' in self.contents[h]:
                g = h + 1
                while g <= len(self.contents)-1:
                    line = self.contents[g].split()
                    if len(line) == 8:
                        relative_energy.append(float(line[1]))
                        total_energy.append(float(line[2]))
                        population.append(float(line[4]))
                        ready = True

                    if '/K' in line[1]:
                        break

                    g += 1
            if ready == True:
                break

            h += 1

        return {'relative energies': relative_energy,
                'total energies': total_energy,
                'population': population}

    def _crest_timing(self):

        ready = False
        for line in self.contents:
            if "test MD wall time" in line:
                test_MD = line
                ready = True

            if "MTD wall time" in line:
                MTD = line

            if "multilevel OPT wall time" in line:
                multilevel_OPT = line

            if "MD wall time" in line and ready == True:
                MD = line
                ready = False

            if "GC wall time" in line:
                GC = line

            if "Overall wall time" in line:
                overall = line

        return {'test MD wall time': test_MD,
                'metadynamics wall time': MTD_time,
                'multilevel opt wall time': multilevel_OPT,
                'molecular dynamics wall time': MD,
                'genetic z-matrix crossing wall time': GC,
                'overall wall time': overall}

    def _isomer_energy(self):
        complete = False
        relative_energies = []
        total_energies = []
        g = len(self.contents)-1
        while g >= 0:
            if 'structure    ΔE(kcal/mol)   Etot(Eh)' in self.contents[g]:
                h = g + 1
                while h <= len(self.contents)-1:
                    if self.contents[h] != ' \n':
                        line = self.contents[h].split()
                        relative_energies.append(float(line[1]))
                        total_energies.append(float(line[2]))
                    else:
                        complete = True
                        break

            if complete == True:
                break

            g -= 1

        return {'relative energy': relative_energies,
                'total energy': total_energies}

    def _isomer_timing(self):

        def time(LINE):
            line = LINE.split(':')
            hr = line[1].split()
            mn = line[2].split()
            sc = line[3].split()
            hr = (LINE.split(':'))[1].split()
            mn = (LINE.split(':'))[2].split()
            sc = (LINE.split(':'))[3].split()
            return hr + mn + sc

        for line in self.contents:
            if "LMO calc. wall time" in line:
                LMO_time = time(line)

            if "multilevel OPT wall time" in line:
                OPT_time = time(line)

            if "Overall wall time" in line:
                OVERALL_time = time(line)

        return {'local molecular orbital wall time': LMO_time,
                'multilevel opt wall time': OPT_time,
                'overall wall time': OVERALL_time}

    def _opt_energy(self):
        for line in self.contents:
            if 'TOTAL ENERGY' in line:
                energy = line.split()[3] + ' Hartrees'

        return None, energy

    def _opt_timing(self):

        def time(LINE):
            line = LINE.split(',')
            hr = line[1].split()
            mn = line[2].split()
            sc = line[3].split()
            return hr + mn + sc

        tot = False
        scf = False
        anc = False

        for line in self.contents:
            if "wall-time" in line and tot is False:
                total_time = time(line)
                tot = True

            elif "wall-time" in line and scf is False:
                scf_time = time(line)
                scf = True

            if "wall-time" in line and anc is False:
                anc_time = time(line)
                anc = True

        return {'Total wall time': total_time,
                'SCF wall time': scf_time,
                'ANC optimizer wall time': anc_time}

    def _parse_energy(self):

        if self.parse_crest == True:
            return self._crest_energy()
        if self.parse_opt == True:
            return self._opt_energy()
        if self.parse_isomer == True:
            return self._isomer_energy()

    def _parse_timing(self):
        if self.parse_crest == True:
            return self._crest_timing()
        if self.parse_opt == True:
            return self._opt_timing()
        if self.parse_isomer == True:
            return self._isomer_timing()

    def _parse_protocol(self):

        protocol = None

        for line in self.contents:
            if ">" in line:
                protocol = line
            elif "program call" in line:
                protocol = line

        return protocol

    def _separate_xyz(self, FILE):
        '''
        Split .xyz into separate XYZGeometry instances
        '''

        if len(list(pybel.readfile('xyz', FILE))) > 1:
            geom_list = []
            count = 1
            XYZ = FILE.split(".")[0]

            for geom in pybel.readfile('xyz', FILE):
                geom.write("xyz", "%s_%d.xyz" % (XYZ, count))
                geom_list.append("%s_%d.xyz" % (XYZ, count))
                count += 1

            x = [isicle.geom.load(i) for i in geom_list]

        else:
            x = [isicle.geom.load(FILE)]

        return isicle.conformers.ConformationalEnsemble(x)

    def _parse_xyz(self):

        FILE = self.xyz_path
        geom = self._separate_xyz(FILE)

        return geom

    def parse(self, to_parse=['energy']):
        '''Extract relevant information from data'''

        # Check that the file is valid first
        if len(self.contents) == 0:
            raise RuntimeError('No contents to parse: {}'.format(self.path))

        self.parse_crest = False
        self.parse_opt = False
        self.parse_isomer = False

        # Initialize result object to store info
        result = XTBResult()

        if self.path.endswith('xyz'):

            if 'geometry' in to_parse:
                try:
                    self.xyz_path = self.path
                    geom = self._parse_xyz()
                    result.set_geometry(geom)
                except IndexError:
                    pass

        if self.path.endswith('out') or self.path.endswith('log'):
            protocol = self._parse_protocol()
            result.set_protocol(protocol)  # Stored as string

            if 'geometry' in to_parse:
                XYZ = None
                if 'xtb' in protocol:
                    self.parse_opt = True
                    XYZ = 'xtbopt.xyz'
                if 'deprotonate' in protocol:
                    self.parse_isomer = True
                    XYZ = 'deprotonate.xyz'
                elif 'protonate' in protocol:
                    self.parse_isomer = True
                    XYZ = 'protonated.xyz'
                elif 'tautomer' in protocol:
                    self.parse_isomer = True
                    XYZ = 'tautomers.xyz'
                elif 'crest' in protocol:
                    self.parse_crest = True
                    XYZ = 'crest_conformers.xyz'

                if XYZ is None:
                    raise RuntimeError('XYZ file associated with XTB job not available, \
                                        please parse separately.')

                else:
                    temp_dir = os.path.dirname(self.path)
                    self.xyz_path = os.path.join(temp_dir, XYZ)
                    try:
                        geom = self._parse_xyz()
                        result.set_geometry(geom)
                    except IndexError:
                        pass

            if 'timing' in to_parse:

                try:
                    timing = self._parse_timing()
                    result.set_timing(timing)  # Stored as dictionary
                except IndexError:
                    pass

            if 'energy' in to_parse:

                try:
                    energy = self._parse_energy()
                    result.set_energy(energy)  # Stored as dictionary
                except IndexError:
                    pass

        self.result = result
        return result

    def save(self, path: str):
        '''Write parsed object to file'''
        self.result.save(path)
        return
