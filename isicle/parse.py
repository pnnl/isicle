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
    """
    Extract text from an NWChem simulation output file.
    """

    def __init__(self):
        self.contents = None
        self.result = None
        self.path = None

    def load(self, path: str):
        """
        Load in the data file
        """
        with open(path, "r") as f:
            self.contents = f.readlines()
        self.path = path
        return self.contents

    def _parse_geometry(self):
        """
        Add docstring
        """
        search = os.path.dirname(self.path)
        geoms = sorted(glob.glob(os.path.join(search, "*.xyz")))

        if len(geoms) > 0:
            return isicle.io.load(geoms[-1])

        raise Exception

    def _parse_energy(self):
        """
        Add docstring
        """
        # TO DO: Add Initial energy and final energy if different

        # Init
        energy = None

        # Cycle through file
        for line in self.contents:
            if "Total DFT energy" in line:
                # Overwrite last saved energy
                energy = float(line.split()[-1])

        return energy

    def _parse_shielding(self):
        """
        Add docstring
        """
        # Init
        ready = False
        shield_idxs = []
        shield_atoms = []
        shields = []

        for line in self.contents:
            if " SHIELDING" in line:
                shield_idxs = [int(x) for x in line.split()[2:]]
                if len(shield_idxs) == 0:
                    collect_idx = True

            if "Atom:" in line:
                atom = line.split()[2]
                idx = line.split()[1]
                ready = True

            elif "isotropic" in line and ready is True:
                shield = float(line.split()[-1])
                shield_atoms.append(atom)
                shields.append(shield)
                if collect_idx is True:
                    shield_idxs.append(int(idx))

        if len(shields) > 1:
            return {"index": shield_idxs, "atom": shield_atoms, "shielding": shields}

        raise Exception

    def _parse_spin(self):
        """
        Add docstring
        """
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

        return {
            "pair indices": coup_pairs,
            "spin couplings": coup,
            "index": index,
            "g-tensors": g_factor,
        }

    def _parse_frequency(self):
        """
        Add docstring
        """
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
            if ("geometry" in line) and (natoms is None):
                atom_start = i + 7
            if ("Atomic Mass" in line) and (natoms is None):
                atom_stop = i - 2
                natoms = atom_stop - atom_start + 1
            if "Normal Eigenvalue" in line:
                has_frequency = True
                freq_start = i + 3
                freq_stop = i + 2 + 3 * natoms

            # Get values
            if "Zero-Point correction to Energy" in line:
                zpe = line.rstrip().split("=")[-1]

            if "Thermal correction to Enthalpy" in line:
                enthalpies = line.rstrip().split("=")[-1]

            if "Total Entropy" in line:
                entropies = line.rstrip().split("=")[-1]

            if "constant volume heat capacity" in line:
                capacities = line.rstrip().split("=    ")[-1]

        if has_frequency is True:
            freq = np.array(
                [float(x.split()[1]) for x in self.contents[freq_start : freq_stop + 1]]
            )
            intensity_au = np.array(
                [float(x.split()[3]) for x in self.contents[freq_start : freq_stop + 1]]
            )
            intensity_debyeangs = np.array(
                [float(x.split()[4]) for x in self.contents[freq_start : freq_stop + 1]]
            )
            intensity_KMmol = np.array(
                [float(x.split()[5]) for x in self.contents[freq_start : freq_stop + 1]]
            )
            intensity_arbitrary = np.array(
                [float(x.split()[6]) for x in self.contents[freq_start : freq_stop + 1]]
            )

            return {
                "frequencies": freq,
                "intensity atomic units": intensity_au,
                "intensity (debye/angs)**2": intensity_debyeangs,
                "intensity (KM/mol)": intensity_KMmol,
                "intensity arbitrary": intensity_arbitrary,
                "correction to enthalpy": enthalpies,
                "total entropy": entropies,
                "constant volume heat capacity": capacities,
                "zero-point correction": zpe,
            }

        raise Exception

    def _parse_charge(self):
        """
        Add docstring
        """
        # TO DO: Parse molecular charge and atomic charges
        # TO DO: Add type of charge
        # TO DO: Multiple instances of charge analysis seen (two Mulliken and one Lowdin, difference?)
        charges = []
        ready = False

        for line in self.contents:
            # Load charges from table
            if "Atom       Charge   Shell Charges" in line:
                # Table header found. Overwrite anything saved previously
                ready = True
                charges = []
            elif ready is True and line.strip() in ["", "Line search:"]:
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
            df = pd.DataFrame(
                [x.split()[0:4] for x in charges],
                columns=["idx", "Atom", "Number", "Charge"],
            )
            df.Number = df.Number.astype("int")
            df.Charge = df.Number - df.Charge.astype("float")

            return df.Charge.values

        raise Exception

    def _parse_timing(self):
        """
        Add docstring
        """
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
            if "No." in line and len(indices) == 0:
                indices.append(i + 2)  # 0
            elif "Atomic Mass" in line and len(indices) == 1:
                indices.append(i - 1)  # 1
                indices.append(i + 3)  # 2
            elif "Effective nuclear repulsion energy" in line and len(indices) == 3:
                indices.append(i - 2)  # 3

            # Check for optimization and frequency calcs
            if "NWChem geometry Optimization" in line:
                opt = True
            elif "NWChem Nuclear Hessian and Frequency Analysis" in line:
                freq = True

            # Get timing
            if "Total iterative time" in line and opt is False:
                preoptTime += float(line.rstrip().split("=")[1].split("s")[0])
            elif "Total iterative time" in line and opt is True and freq is False:
                geomoptTime += float(line.rstrip().split("=")[1].split("s")[0])
            elif "Total iterative time" in line and freq is True:
                freqTime += float(line.rstrip().split("=")[1].split("s")[0])

            if "Total times" in line:
                cpuTime = float(line.rstrip().split(":")[1].split("s")[0])
                # wallTime = float(line.rstrip().split(':')[2].split('s')[0])
                freqTime = cpuTime - geomoptTime - preoptTime

        # natoms = int(self.contents[indices[1] - 1].split()[0])

        return {
            "single point": preoptTime,
            "geometry optimization": geomoptTime,
            "frequency": freqTime,
            "total": cpuTime,
        }

    def _parse_molden(self):
        """
        Add docstring
        """
        search = splitext(self.path)[0]
        m = glob.glob(search + "*.molden")

        if len(m) != 1:
            raise Exception

        return m[0]

    def _parse_protocol(self):
        """
        Parse out dft protocol
        """
        functional = []
        basis_set = []
        solvation = []
        tasks = []
        basis = None
        func = None
        solvent = None

        for line in self.contents:
            if "* library" in line:
                basis = line.split()[-1]
            if " xc " in line:
                func = line.split(" xc ")[-1].strip()
            if "solvent " in line:
                solvent = line.split()[-1]
            if "task dft optimize" in line:
                tasks.append("optimize")
                basis_set.append(basis)
                functional.append(func)
                solvation.append(solvent)
            if "SHIELDING" in line:
                tasks.append("shielding")
                basis_set.append(basis)
                functional.append(func)
                solvation.append(solvent)
            if "SPINSPIN" in line:
                tasks.append("spin")
                basis_set.append(basis)
                functional.append(func)
                solvation.append(solvent)
            if "freq " in line:
                tasks.append("frequency")
                basis_set.append(basis)
                functional.append(func)
                solvation.append(solvent)

        return {
            "functional": functional,
            "basis set": basis_set,
            "solvation": solvation,
            "tasks": tasks,
        }

    def _parse_connectivity(self):
        """
        Add docstring
        """
        coor_substr = "internuclear distances"

        # Extracting Atoms & Coordinates
        ii = [i for i in range(len(self.contents)) if coor_substr in self.contents[i]]
        ii.sort()

        g = ii[0] + 4
        connectivity = []
        while g <= len(self.contents) - 1:
            if "-" not in self.contents[g]:
                line = self.contents[g].split()
                pair = [line[1], line[4], int(line[0]), int(line[3])]
                connectivity.append(pair)

            else:
                break
            g += 1

        return connectivity

    def parse(self):
        """
        Extract relevant information from NWChem output

        Parameters
        ----------
        to_parse : list of str
            geometry, energy, shielding, spin, frequency, molden, charge, timing
        """

        # Check that the file is valid first
        if len(self.contents) == 0:
            raise RuntimeError("No contents to parse: {}".format(self.path))
        if "Total times  cpu" not in self.contents[-1]:
            raise RuntimeError("Incomplete NWChem run: {}".format(self.path))

        # Initialize result object to store info
        result = {}

        try:
            result["protocol"] = self._parse_protocol()
        except:
            pass

        try:
            result["geom"] = self._parse_geometry()

        except:
            pass

        try:
            result["energy"] = self._parse_energy()
        except:
            pass

        try:
            result["shielding"] = self._parse_shielding()
        except:  # Must be no shielding info
            pass

        try:
            result["spin"] = self._parse_spin()
        except:
            pass

        try:
            result["frequency"] = self._parse_frequency()
        except:
            pass

        try:
            result["molden"] = self._parse_molden()
        except:
            pass

        try:
            result["charge"] = self._parse_charge()
        except:
            pass

        try:
            result["timing"] = self._parse_timing()
        except:
            pass

        try:
            result["connectivity"] = self._parse_connectivity()
        except:
            pass

        return result

    def save(self, path):
        with open(path, "wb") as f:
            pickle.dump(self, f)
        return


class ImpactParser(FileParserInterface):
    """
    Extract text from an Impact mobility calculation output file.
    """

    def __init__(self):
        """
        Add docstring
        """
        self.contents = None
        self.result = None

    def load(self, path: str):
        """
        Load in the data file
        """
        with open(path, "rb") as f:
            self.contents = f.readlines()

        return self.contents

    def parse(self):
        """
        Extract relevant information from data
        """

        # Check CCS results == 1
        count = 0
        for line in self.contents:
            l = line.split(" ")
            if "CCS" in l[0]:
                count += 1
        if count != 1:
            return self.result

        # Assume values in second line
        l = self.contents[1].split(" ")
        l = [x for x in l if len(x) > 0]

        # Pull values of interest - may be error prone
        values = []
        try:
            values.append(float(l[-5]))
            values.append(float(l[-3][:-1]))
            values.append(float(l[-2]))
            values.append(int(l[-1]))
        except (ValueError, IndexError) as e:
            print("Could not parse file: ", e)
            return None

        # Add to dictionary to return
        result = {}
        keys = ["CCS_PA", "SEM_rel", "CCS_TJM", "n_iter"]
        for key, val in zip(keys, values):
            result[key] = [val]

        # Save and return results
        self.result = result
        return result  # TODO: return CCS?

    def save(self, path: str, sep="\t"):
        """
        Write parsed object to file
        """
        pd.DataFrame(self.result).to_csv(path, sep=sep, index=False)
        return


class MobcalParser(FileParserInterface):
    """
    Extract text from a MOBCAL mobility calculation output file.
    """

    def __init__(self):
        """
        Add docstring
        """
        self.contents = None
        self.result = {}

    def load(self, path: str):
        """
        Load in the data file
        """
        with open(path, "r") as f:
            self.contents = f.readlines()

        return self.contents

    def parse(self):
        """
        Extract relevant information from data
        """
        done = False
        for line in self.contents:
            # if "average (second order) TM mobility" in line:
            #     m_mn = float(line.split('=')[-1])
            if "average TM cross section" in line:
                ccs_mn = float(line.split("=")[-1])
            elif "standard deviation TM cross section" in line:
                ccs_std = float(line.split("=")[-1])
            elif "standard deviation (percent)" in line:
                done = True
        if done is True:
            self.result["ccs"] = {"mean": ccs_mn, "std": ccs_std}

        return self.result

    def save(self, path: str, sep="\t"):
        """
        Write parsed object to file
        """
        pd.DataFrame(self.result).to_csv(path, sep=sep, index=False)
        return


class SanderParser(FileParserInterface):
    """
    Extract text from an Sander simulated annealing simulation output file.
    """

    def load(self, path: str):
        """
        Load in the data file
        """
        raise NotImplementedError

    def parse(self):
        """
        Extract relevant information from data
        """
        raise NotImplementedError

    def save(self, path: str):
        """
        Write parsed object to file
        """
        raise NotImplementedError


class XTBParser(FileParserInterface):
    """
    Add docstring
    """

    def __init__(self):
        """
        Add docstring
        """
        self.contents = None
        self.result = None
        self.path = None

    def load(self, path: str):
        """
        Load in the data file
        """
        with open(path, "r") as f:
            self.contents = f.readlines()
        self.path = path
        # return self.contents

    def _crest_energy(self):
        """
        Add docstring
        """
        relative_energy = []
        total_energy = []
        population = []

        ready = False
        for h in range(len(self.contents), 0, -1):
            if "Erel/kcal" in self.contents[h]:
                g = h + 1
                for j in range(g, len(self.contents)):
                    line = self.contents[j].split()
                    if len(line) == 8:
                        relative_energy.append(float(line[1]))
                        total_energy.append(float(line[2]))
                        population.append(float(line[4]))
                        ready = True

                    if "/K" in line[1]:
                        break
            if ready == True:
                break

        return {
            "relative energies": relative_energy,
            "total energies": total_energy,
            "population": population,
        }

    def _crest_timing(self):
        """
        Add docstring
        """

        def grab_time(line):
            line = line.replace(" ", "")
            line = line.split(":")

            return ":".join(line[1:]).strip("\n")

        ready = False
        for line in self.contents:
            if "test MD wall time" in line:
                test_MD = grab_time(line)
                ready = True

            if "MTD wall time" in line:
                MTD_time = grab_time(line)

            if "multilevel OPT wall time" in line:
                multilevel_OPT = grab_time(line)

            if "MD wall time" in line and ready == True:
                MD = grab_time(line)
                ready = False

            if "GC wall time" in line:
                GC = grab_time(line)

            if "Overall wall time" in line:
                overall = grab_time(line)

        return {
            "test MD wall time": test_MD,
            "metadynamics wall time": MTD_time,
            "multilevel opt wall time": multilevel_OPT,
            "molecular dynamics wall time": MD,
            "genetic z-matrix crossing wall time": GC,
            "overall wall time": overall,
        }

    def _isomer_energy(self):
        """
        Add docstring
        """
        complete = False
        relative_energies = []
        total_energies = []
        for i in range(len(self.contents), 0, -1):
            if "structure    ΔE(kcal/mol)   Etot(Eh)" in self.contents[i]:
                h = i + 1
                for j in range(h, len(self.contents)):
                    if self.contents[j] != " \n":
                        line = self.contents[j].split()
                        relative_energies.append(float(line[1]))
                        total_energies.append(float(line[2]))
                    else:
                        complete = True
                        break

            if complete == True:
                break

        return {"relative energy": relative_energies, "total energy": total_energies}

    def _isomer_timing(self):
        """
        Add docstring
        """

        def grab_time(line):
            line = line.replace(" ", "")
            line = line.split(":")

            return ":".join(line[1:]).strip("\n")

        for line in self.contents:
            if "LMO calc. wall time" in line:
                LMO_time = grab_time(line)

            if "multilevel OPT wall time" in line:
                OPT_time = grab_time(line)

            if "Overall wall time" in line:
                OVERALL_time = grab_time(line)

        return {
            "local molecular orbital wall time": LMO_time,
            "multilevel opt wall time": OPT_time,
            "overall wall time": OVERALL_time,
        }

    def _opt_energy(self):
        """
        Add docstring
        """
        for line in self.contents:
            if "TOTAL ENERGY" in line:
                energy = line.split()[3] + " Hartrees"

        return {"Total energy": energy}

    def _opt_timing(self):
        """
        Add docstring
        """

        def grab_time(line):
            line = line.replace(" ", "")
            line = line.split(":")

            return ":".join(line[1:]).strip("\n")

        tot = False
        scf = False
        anc = False

        for line in self.contents:
            if "wall-time" in line and tot is False:
                total_time = grab_time(line)
                tot = True

            elif "wall-time" in line and scf is False:
                scf_time = grab_time(line)
                scf = True

            if "wall-time" in line and anc is False:
                anc_time = grab_time(line)
                anc = True

        return {
            "Total wall time": total_time,
            "SCF wall time": scf_time,
            "ANC optimizer wall time": anc_time,
        }

    def _parse_energy(self):
        """
        Add docstring
        """
        if self.parse_crest == True:
            return self._crest_energy()
        if self.parse_opt == True:
            return self._opt_energy()
        if self.parse_isomer == True:
            return self._isomer_energy()

    def _parse_timing(self):
        """
        Add docstring
        """
        if self.parse_crest == True:
            return self._crest_timing()
        if self.parse_opt == True:
            return self._opt_timing()
        if self.parse_isomer == True:
            return self._isomer_timing()

    def _parse_protocol(self):
        """
        Add docstring
        """
        protocol = None

        for line in self.contents:
            if " > " in line:
                protocol = line.strip("\n")
            if "program call" in line:
                protocol = (line.split(":")[1]).strip("\n")
        return protocol

    def _parse_xyz(self):
        """
        Split .xyz into separate XYZGeometry instances
        """

        FILE = self.xyz_path
        if len(list(pybel.readfile("xyz", FILE))) > 1:
            geom_list = []
            count = 1
            XYZ = FILE.split(".")[0]
            for geom in pybel.readfile("xyz", FILE):
                geom.write("xyz", "%s_%d.xyz" % (XYZ, count))
                geom_list.append("%s_%d.xyz" % (XYZ, count))
                count += 1

            x = [isicle.io.load(i) for i in geom_list]

        else:
            x = [isicle.io.load(self.xyz_path)]

        return isicle.conformers.ConformationalEnsemble(x)

    def parse(self):
        """
        Extract relevant information from data
        """

        # Check that the file is valid first
        if len(self.contents) == 0:
            raise RuntimeError("No contents to parse: {}".format(self.path))
        if "terminated normally" not in self.contents[-1]:
            if "ratio" not in self.contents[-2]:
                raise RuntimeError("XTB job failed: {}".format(self.path))

        self.parse_crest = False
        self.parse_opt = False
        self.parse_isomer = False

        # Initialize result object to store info
        result = {}

        try:
            result["protocol"] = self._parse_protocol()
        except:
            pass
        try:
            result["timing"] = self._parse_timing()
        except:
            pass

        try:
            result["energy"] = self._parse_energy()
        except:
            pass

        # Parse geometry from assoc. XYZ file
        try:
            if self.path.endswith("xyz"):
                if "geometry" in to_parse:
                    try:
                        self.xyz_path = self.path
                        result["geom"] = self._parse_xyz()

                    except:
                        pass

            if self.path.endswith("out") or self.path.endswith("log"):
                # try geometry parsing
                try:
                    XYZ = None
                    if result["protocol"].split()[0] == "xtb":
                        self.parse_opt = True
                        XYZ = "xtbopt.xyz"
                    if result["protocol"].split()[1] == "crest":
                        if "-deprotonate" in result["protocol"]:
                            self.parse_isomer = True
                            XYZ = "deprotonated.xyz"
                        elif "-protonate" in result["protocol"]:
                            self.parse_isomer = True
                            XYZ = "protonated.xyz"
                        elif "-tautomer" in result["protocol"]:
                            self.parse_isomer = True
                            XYZ = "tautomers.xyz"
                        else:
                            self.parse_crest = True
                            XYZ = "crest_conformers.xyz"

                    if XYZ is None:
                        raise RuntimeError(
                            "XYZ file associated with XTB job not available,\
                                        please parse separately."
                        )

                    else:
                        temp_dir = os.path.dirname(self.path)
                        self.xyz_path = os.path.join(temp_dir, XYZ)

                        result["geom"] = self._parse_xyz()
                except:
                    pass
        except:
            pass
        return result

    def save(self, path):
        """
        Add docstring
        """
        with open(path, "wb") as f:
            pickle.dump(self, f)
        return


class RDKitParser(FileParserInterface):
    """ """

    def __init__(self):
        """
        Add docstring
        """
        return

    def load(self, path: str):
        """
        Add docstring
        """
        return

    def parse(self):
        """
        Add docstring
        """
        return

    def save(self, path: str):
        """
        Add docstring
        """
        return


class TINKERParser(FileParserInterface):
    """
    Add docstring
    """

    def __init__(self):
        """
        Add docstring
        """
        self.contents = None
        self.result = None
        self.path = None

    def load(self, path: str):
        """
        Load in the data file
        """
        with open(path, "r") as f:
            self.contents = f.readlines()
        self.path = path

    def _parse_energy(self):
        """
        Add docstring
        """
        inp = self.contents
        if len(inp) < 13:
            quit()

        # Get the conformer energies from the file
        energies = []
        for line in inp[13:]:
            data = line[:-1].split("  ")
            data = [_f for _f in data if _f]
            if len(data) >= 3:
                if "Map" in data[0] and "Minimum" in data[1]:
                    energies.append(float(data[-1]))

        return energies

    def _parse_conformers(self):
        """
        Add docstring
        """

        def parse_atom_symbol(AtomNum):
            Lookup = [
                "H",
                "He",
                "Li",
                "Be",
                "B",
                "C",
                "N",
                "O",
                "F",
                "Ne",
                "Na",
                "Mg",
                "Al",
                "Si",
                "P",
                "S",
                "Cl",
                "Ar",
                "K",
                "Ca",
                "Sc",
                "Ti",
                "V",
                "Cr",
                "Mn",
                "Fe",
                "Co",
                "Ni",
                "Cu",
                "Zn",
                "Ga",
                "Ge",
                "As",
                "Se",
                "Br",
                "Kr",
                "Rb",
                "Sr",
                "Y",
                "Zr",
                "Nb",
                "Mo",
                "Tc",
                "Ru",
                "Rh",
                "Pd",
                "Ag",
                "Cd",
                "In",
                "Sn",
                "Sb",
                "Te",
                "I",
                "Xe",
                "Cs",
                "Ba",
                "La",
                "Ce",
                "Pr",
                "Nd",
                "Pm",
                "Sm",
                "Eu",
                "Gd",
                "Tb",
                "Dy",
                "Ho",
                "Er",
                "Tm",
                "Yb",
                "Lu",
                "Hf",
                "Ta",
                "W",
                "Re",
                "Os",
                "Ir",
                "Pt",
                "Au",
                "Hg",
                "Tl",
                "Pb",
                "Bi",
                "Po",
                "At",
                "Rn",
            ]

            if AtomNum > 0 and AtomNum < len(Lookup):
                return Lookup[AtomNum - 1]
            else:
                print("No such element with atomic number " + str(AtomNum))
                return 0

        conffile = open(self.path.split(".")[0] + ".arc", "r")
        confdata = conffile.readlines()
        conffile.close()
        conformers = []
        atoms = []
        atomtypes = isicle.utils.tinker_lookup()["atomtypes"].to_list()
        anums = isicle.utils.tinker_lookup()["anums"].to_list()
        atypes = [x[:3] for x in atomtypes]

        # Parse data from arc file, extract coordinates and atom type
        for line in confdata:
            data = [_f for _f in line.split("  ") if _f]
            if len(data) < 3:
                conformers.append([])
            else:
                if len(conformers) == 1:
                    anum = anums[atypes.index(data[1][:3])]
                    atoms.append(parse_atom_symbol(anum))
                conformers[-1].append([x for x in data[2:5]])

        # Convert from TINKER xyz format to standard xyz format
        xyz_file = []
        for conf in conformers:
            xyz_file.append(" {}\n".format(len(conf)))
            for idx, line in enumerate(conf):
                s = " {}\t{}\t{}\t{}".format(atoms[idx], line[0], line[1], line[2])
                xyz_file.append(s)

        # Write xyzs to file
        FILE = "conformers.xyz"
        f = open(FILE, "w+")
        for i in xyz_file:
            f.write(i + "\n")
        f.close()

        # Read in file by
        if len(list(pybel.readfile("xyz", FILE))) > 1:
            geom_list = []
            count = 1
            XYZ = FILE.split(".")[0]
            for geom in pybel.readfile("xyz", FILE):
                geom.write("xyz", "%s_%d.xyz" % (XYZ, count), overwrite=True)
                geom_list.append("%s_%d.xyz" % (XYZ, count))
                count += 1

            x = [isicle.io.load(i) for i in geom_list]

        else:
            x = [isicle.io.load(self.xyz_path)]

        return isicle.conformers.ConformationalEnsemble(x)

    def parse(self):
        """
        Extract relevant information from data
        """

        # Check that the file is valid first
        if len(self.contents) == 0:
            raise RuntimeError("No contents to parse: {}".format(self.path))

        # Initialize result object to store info
        result = {}

        try:
            result["geom"] = self._parse_conformers()
        except:
            pass

        try:
            result["energy"] = self._parse_energy()
        except:
            pass

        try:
            result["charge"] = self._parse_charge()
        except:
            pass

        return result

    def save(self):
        """
        Add docstring
        """
        return
