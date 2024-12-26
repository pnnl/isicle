import pickle
import re

import numpy as np
import pandas as pd

import isicle
from isicle.interfaces import FileParserInterface


class ORCAParser(FileParserInterface):
    """Extract information from an ORCA simulation output files."""

    def __init__(self, data=None):
        self.data = data

        self.result = {}

    def load(self, path):
        self.data = isicle.io.load_pickle(path)

    def _find_output_by_header(self, header):
        # Fat regex
        pattern = (
            r"(-{2,})\n\s{0,}("
            + re.escape(header)
            + r")\s{0,}\n-{2,}\n([\s\S]*?)(?=-{2,}\n[^\n]*\n-{2,}\n|$)"
        )

        # Search ORCA output file
        matches = re.findall(pattern, self.data["out"])

        # Return "body" of each match
        return [x[2].strip() for x in matches]

    def _parse_protocol(self):
        return self.data["inp"]

    def _parse_geometry(self):
        return self.data["xyz"]

    def _parse_energy(self):
        # Split text
        lines = self.data["property"].split("\n")

        # Search for energy values
        elines = [x for x in lines if "Total DFT Energy" in x]

        # Energy values not found
        if len(elines) == 0:
            return None

        # Map float over values
        evals = [float(x.split()[-1].strip()) for x in elines]

        # Return last energy value
        return evals[-1]

    def _parse_frequency(self):
        if "hess" in self.data:
            # Define columns
            columns = ["wavenumber", "eps", "intensity", "TX", "TY", "TZ"]

            # Split sections by delimiter
            blocks = self.data["hess"].split("$")

            # Search for frequency values
            freq_block = [x for x in blocks if x.startswith("ir_spectrum")]

            # Frequency values not found
            if len(freq_block) == 0:
                return None

            # Grab last frequency block
            # Doubtful if more than one, but previous results in list
            freq_block = freq_block[-1]

            # Split block into lines
            lines = freq_block.split("\n")

            # Map float over values
            vals = np.array(
                [
                    list(map(float, x.split()))
                    for x in lines
                    if len(x.split()) == len(columns)
                ]
            )

            # Zip columns and values
            return dict(zip(columns, vals.T))
        
        # No frequency info
        return None

    def _parse_timing(self):
        # Grab only last few lines
        lines = self.data["out"].split("\n")[-100:]

        # Find start of timing section
        parts = []
        start_idx = None
        for i, line in enumerate(lines):
            if line.startswith("Timings for individual modules"):
                start_idx = i + 2

            # Strip out extraneous info
            parts.append(
                [x.strip() for x in line.split("  ") if x and x.strip() != "..."]
            )

        # Timing not found
        if start_idx is None:
            return None

        # Split out timing section
        tlines = lines[start_idx:]
        tparts = parts[start_idx:]

        # Individual timings
        timings = [x for x in tparts if any([" sec " in y for y in x])]
        timings = {x[0].strip("..."): float(x[1].split()[0]) for x in timings}

        # Boolean indication of success
        success = len([x for x in tlines if "ORCA TERMINATED NORMALLY" in x]) > 0
        timings["success"] = success

        # Total time
        total_time = [x for x in tlines if "TOTAL RUN TIME" in x]

        if len(total_time) > 0:
            total_time = total_time[-1].split(":")[-1].strip()
            times = list(map(int, total_time.split()[::2]))
            units = total_time.split()[1::2]
        else:
            total_time = None

        timings["Total run time"] = dict(zip(units, times))

        return timings

    def _parse_shielding(self):
        # Filter comments
        property = [
            x.strip()
            for x in self.data["property"].split("\n")
            if not x.startswith("#")
        ]
        property = "\n".join(property)

        # Split sections by delimiter
        blocks = property.split("$ ")

        # Search for shielding values
        shielding_block = [x for x in blocks if x.startswith("EPRNMR_OrbitalShielding")]

        # Shielding values not found
        if len(shielding_block) == 0:
            return None

        # Grab last shielding block
        # Doubtful if more than one, but previous results in list
        shielding_block = shielding_block[-1]

        # Define a pattern for extracting relevant information
        pattern = re.compile(
            r"Nucleus: (\d+) (\w+)\n(Shielding tensor.*?P\(iso\) \s*[-+]?\d*\.\d+)",
            re.DOTALL,
        )

        # Match against pattern
        matches = pattern.findall(shielding_block)

        # Result container
        shielding = {}

        # Enumerate matches
        for match in matches:
            # Per-nucleus info
            nucleus_index = match[0]
            nucleus_name = match[1]
            nucleus_data = match[2]

            # Extracting values using regex
            tensors = re.findall(r"(-?\d+\.\d+|-?\d+.\d+e[+-]\d+)", nucleus_data)
            tensors = [float(val) for val in tensors]

            # Creating arrays from extracted values
            shielding_tensor = np.array(tensors[:9]).reshape(3, 3)
            p_tensor_eigenvectors = np.array(tensors[9:18]).reshape(3, 3)
            p_eigenvalues = np.array(tensors[18:21])
            p_iso = float(tensors[21])

            # Constructing the dictionary with nuclei index and name
            shielding[f"{nucleus_index}{nucleus_name}"] = {
                "shielding tensor": shielding_tensor,
                "P tensor eigenvectors": p_tensor_eigenvectors,
                "P eigenvalues": p_eigenvalues,
                "P(iso)": p_iso,
            }

        # Add shielding summary
        shielding["shielding_summary"] = self._parse_shielding_summary()

        return shielding

    def _parse_orbital_energies(self):
        header = "ORBITAL ENERGIES"
        text = self._find_output_by_header(header)

        # Orbital energies not found
        if len(text) == 0:
            return None

        # Get last relevant output
        text = text[-1].split("\n")

        # Parse table
        text = [x.strip() for x in text if x.strip() != "" and "*" not in x]
        columns = text[0].split()
        body = [x.split() for x in text[1:]]

        # Construct data frame
        df = pd.DataFrame(body, columns=columns, dtype=float)

        # Map correct types
        df["NO"] = df["NO"].astype(int)

        # Drop unoccupied orbitals?
        return df

    def _parse_spin(self):
        header = "SUMMARY OF ISOTROPIC COUPLING CONSTANTS (Hz)"
        text = self._find_output_by_header(header)

        # Spin couplings not found
        if len(text) == 0:
            return None

        # Get last relevant output
        text = text[-1].split("\n")

        # Parse table
        text = [x.strip() for x in text if x.strip() != "" and "*" not in x]
        columns = [x.replace(" ", "") for x in re.split(r"\s{2,}", text[0])]
        body = [re.split(r"\s{2,}", x)[1:] for x in text[1:-1]]

        # Construct data frame
        return pd.DataFrame(body, dtype=float, columns=columns, index=columns)

    def _parse_shielding_summary(self):
        header = "CHEMICAL SHIELDING SUMMARY (ppm)"
        text = self._find_output_by_header(header)

        # Shielding values not found
        if len(text) == 0:
            return None

        # Get last relevant output
        text = text[-1].split("\n")

        # Parse table
        text = [x.strip() for x in text if x.strip() != ""]

        # Find stop index
        stop_idx = -1
        for i, row in enumerate(text):
            if all([x == "-" for x in row]):
                stop_idx = i
                break

        # Split columns and body
        columns = text[0].split()
        body = [x.split() for x in text[2:stop_idx]]

        # Construct data frame
        df = pd.DataFrame(body, columns=columns)

        # Map correct types
        for col, dtype in zip(df.columns, (int, str, float, float)):
            df[col] = df[col].astype(dtype)
        return df

    def _parse_thermo(self):
        # In hessian file
        header = "THERMOCHEMISTRY_Energies"

    def _parse_molden(self):
        return None

    def _parse_charge(self):
        return None

    def _parse_connectivity(self):
        return None

    def parse(self):
        result = {
            "protocol": self._parse_protocol(),
            "geometry": self._parse_geometry(),
            "total_dft_energy": self._parse_energy(),
            "orbital_energies": self._parse_orbital_energies(),
            "shielding": self._parse_shielding(),
            "spin": self._parse_spin(),
            "frequency": self._parse_frequency(),
            "molden": self._parse_molden(),
            "charge": self._parse_charge(),
            "timing": self._parse_timing(),
            "connectivity": self._parse_connectivity(),
        }

        # Pop success from timing
        if result["timing"] is not None:
            result["success"] = result["timing"].pop("success")
        else:
            result["success"] = False

        # Filter empty fields
        result = {k: v for k, v in result.items() if v is not None}

        # Add result info to geometry object
        if "geometry" in result:
            result["geometry"].add___dict__(
                {k: v for k, v in result.items() if k != "geometry"}
            )

        # Store attribute
        self.result = result

        return result

    def save(self, path):
        isicle.io.save_pickle(path, self.result)


class NWChemParser(FileParserInterface):
    """Extract information from NWChem simulation output files."""

    def __init__(self, data=None):
        self.data = data

        self.result = {}

    def load(self, path):
        self.data = isicle.io.load_pickle(path)

    def _parse_geometry(self):
        """
        Add docstring
        """

        return self.data["xyz"]["final"]

    def _parse_energy(self):
        """
        Add docstring
        """
        # TO DO: Add Initial energy and final energy if different

        # Init
        energy = None

        # Cycle through file
        for line in self.data["out"].split("\n"):
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

        for line in self.data["out"].split("\n"):
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

        # No shielding data found
        return None

    def _parse_spin(self):
        """
        Add docstring
        """

        # No spin
        if "SPINSPIN" not in self.data["nw"]:
            return None

        # TO DO: Add g-factors

        # Declaring couplings
        coup_pairs = []
        coup = []
        index = []
        g_factor = []
        ready = False

        for line in self.data["out"].split("\n"):
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

        if len(coup_pairs) > 0:
            return {
                "pair indices": coup_pairs,
                "spin couplings": coup,
                "index": index,
                "g-tensors": g_factor,
            }

        # No spin data found
        return None

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

        lines = self.data["out"].split("\n")
        for i, line in enumerate(lines):
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
                [float(x.split()[1]) for x in lines[freq_start : freq_stop + 1]]
            )
            intensity_au = np.array(
                [float(x.split()[3]) for x in lines[freq_start : freq_stop + 1]]
            )
            intensity_debyeangs = np.array(
                [float(x.split()[4]) for x in lines[freq_start : freq_stop + 1]]
            )
            intensity_KMmol = np.array(
                [float(x.split()[5]) for x in lines[freq_start : freq_stop + 1]]
            )
            intensity_arbitrary = np.array(
                [float(x.split()[6]) for x in lines[freq_start : freq_stop + 1]]
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

        # No frequency data found
        return None

    def _parse_charge(self):
        """
        Add docstring
        """
        # TO DO: Parse molecular charge and atomic charges
        # TO DO: Add type of charge
        # TO DO: Multiple instances of charge analysis seen (two Mulliken and one Lowdin, difference?)
        charges = []
        ready = False

        for line in self.data["out"].split("\n"):
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

        # No charge data found
        return None

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

        for i, line in enumerate(self.data["out"].split("\n")):
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

        if cpuTime != 0:
            return {
                "single point": preoptTime,
                "geometry optimization": geomoptTime,
                "frequency": freqTime,
                "total": cpuTime,
                "success": True,
            }
        return None

    def _parse_molden(self):
        """
        Add docstring
        """
        if "molden" in self.data:
            return self["molden"]

        return None

    def _parse_protocol(self):
        """
        Parse out dft protocol
        """
        return self.data["nw"]

    def _parse_connectivity(self):
        """
        Add docstring
        """

        # Split lines
        lines = self.data["out"].split("\n")

        # Extracting Atoms & Coordinates
        coor_substr = "internuclear distances"
        ii = [i for i in range(len(lines)) if coor_substr in lines[i]]

        # Exit condition
        if len(ii) == 0:
            return None

        # Sort hits
        ii.sort()

        # Iterate connectivity
        g = ii[0] + 4
        connectivity = []
        while g <= len(lines) - 1:
            if "-" not in lines[g]:
                line = lines[g].split()
                pair = [line[1], line[4], int(line[0]), int(line[3])]
                connectivity.append(pair)

            else:
                break
            g += 1

        # Check for result
        if len(connectivity) > 0:
            return connectivity

        return None

    def parse(self):
        result = {
            "protocol": self._parse_protocol(),
            "geometry": self._parse_geometry(),
            "total_dft_energy": self._parse_energy(),
            # "orbital_energies": self._parse_orbital_energies(),
            "shielding": self._parse_shielding(),
            "spin": self._parse_spin(),
            "frequency": self._parse_frequency(),
            "molden": self._parse_molden(),
            "charge": self._parse_charge(),
            "timing": self._parse_timing(),
            "connectivity": self._parse_connectivity(),
        }

        # Pop success from timing
        if result["timing"] is not None:
            result["success"] = result["timing"].pop("success")
        else:
            result["success"] = False

        # Filter empty fields
        result = {k: v for k, v in result.items() if v is not None}

        # Add result info to geometry object
        if "geometry" in result:
            result["geometry"].add___dict__(
                {k: v for k, v in result.items() if k != "geometry"}
            )

        # Store attribute
        self.result = result

        return result

    def save(self, path):
        isicle.io.save_pickle(path, self.result)


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

    def __init__(self, data=None):
        self.data = data

        self.result = {}

        if data is not None:
            self.lines = self.data["out"].split("\n")

    def load(self, path):
        self.data = isicle.io.load_pickle(path)

        self.lines = self.data["out"].split("\n")

    def _crest_energy(self):
        """
        Add docstring
        """
        relative_energy = []
        total_energy = []
        population = []

        ready = False
        for h in range(len(self.lines) - 1, -1, -1):
            if "Erel/kcal" in self.lines[h]:
                g = h + 1
                for j in range(g, len(self.lines)):
                    line = self.lines[j].split()
                    if len(line) == 8:
                        relative_energy.append(float(line[1]))
                        total_energy.append(float(line[2]))
                        population.append(float(line[4]))
                        ready = True

                    if "/K" in line[1]:
                        break
            if ready is True:
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
            # Regular expression pattern
            pattern = r"(?:(\d+)\s*d,\s*)?(?:(\d+)\s*h,\s*)?(?:(\d+)\s*min,\s*)?([\d.]+)\s*sec"
            match = re.search(pattern, line)
            return {
                "days": int(match.group(1)) if match.group(1) else 0,
                "hours": int(match.group(2)) if match.group(2) else 0,
                "minutes": int(match.group(3)) if match.group(3) else 0,
                "seconds": float(match.group(4))
                }

        timing = {}
        for line in self.lines:
            if "CREST runtime (total)" in line:
                timing["CREST runtime (total)"] = grab_time(line)

            if "Trial metadynamics (MTD)" in line:
                timing["Trial metadynamics (MTD)"] = grab_time(line)

            if "Metadynamics (MTD)" in line:
                timing["Metadynamics (MTD)"] = grab_time(line)

            if "Geometry optimization" in line:
                timing["Geometry optimization"] = grab_time(line)

            if "Molecular dynamics (MD)" in line:
                timing["Molecular dynamics (MD)"] = grab_time(line)

            if "Genetic crossing (GC)" in line:
                timing["Genetic crossing (GC)"] = grab_time(line)
            
            if "I/O and setup" in line:
                timing["I/O and setup"] = grab_time(line)

        return timing

    def _isomer_energy(self):
        """
        Add docstring
        """
        complete = False
        relative_energies = []
        total_energies = []
        for i in range(len(self.lines) - 1, -1, -1):
            if "structure    ΔE(kcal/mol)   Etot(Eh)" in self.lines[i]:
                h = i + 1
                for j in range(h, len(self.lines)):
                    if self.lines[j] != " \n":
                        line = self.lines[j].split()
                        relative_energies.append(float(line[1]))
                        total_energies.append(float(line[2]))
                    else:
                        complete = True
                        break

            if complete is True:
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

        for line in self.lines:
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
        for line in self.lines:
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

        for line in self.lines:
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

    def _parse_protocol(self):
        """
        Add docstring
        """
        protocol = None

        for line in self.lines:
            if "$ crest" in line:
                protocol = line.strip("\n")
                return protocol
            if "program call" in line:
                protocol = (line.split(":")[1]).strip("\n")
                return protocol
        return protocol

    def _parse_geometry(self):
        """
        Split .xyz into separate XYZGeometry instances
        """
        geometries = {}
        # Add geometry info
        for key in [
            "conformers",
            "rotamers",
            "final",
            "best",
            "protonated",
            "deprotonated",
            "tautomers",
        ]:
            if key in self.data:
                geometries[key] = self.data[key]
        
        if len(geometries) > 1:
            return geometries
        
        return geometries.popitem()[1]

    # TODO
    def _parse_orbital_energies(self):
        pass

    def parse(self):
        """
        Extract relevant information from data
        """

        # Check that the file is valid first
        if len(self.lines) == 0:
            raise RuntimeError("No contents to parse.")

        last_lines = "".join(self.lines[-10:])
        if (
            ("terminat" not in last_lines)
            & ("normal" not in last_lines)
            & ("ratio" not in last_lines)
        ):
            raise RuntimeError("XTB job failed.")

        # Initialize result object to store info
        result = {
            "protocol": self._parse_protocol(),
            "geometry": self._parse_geometry()
        }

        if result["protocol"].split()[0] == "xtb":
            result["timing"] = self._opt_timing()
            result["energy"] = self._opt_energy()

        elif result["protocol"].split()[1] == "crest":
            if any(
                [
                    x in result["protocol"]
                    for x in ["-deprotonate", "-protonate", "-tautomer"]
                ]
            ):
                result["timing"] = self._isomer_timing()
                result["energy"] = self._isomer_energy()
            else:
                result["timing"] = self._crest_timing()
                result["energy"] = self._crest_energy()

        return result

    def save(self, path):
        """
        Add docstring
        """
        with open(path, "wb") as f:
            pickle.dump(self, f)
        return
