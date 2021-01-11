import isicle
from isicle.interfaces import GeometryGenerationInterface
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

## TO DO, read from xyz or mol2 to yield class
## in utils.py, read_mol, Mol, pop_aom, push_atom functions

class GeometryGeneration(GeometryGenerationInterface):
    def __init__(self):
        self.contents = None

    def load(self, path: str):
        """Load in SMILES string"""
        # offer commandline string

        with open(path, 'r') as f:
            self.contents=f.readlines

        m = Chem.MolFromSmiles(self.contents)
        m2 = Chem.AddHs(m)

        return m2
    
    def convert2D(self, path: str):
        """Convert into 2D geometry using RDKit"""

        AllChem.Compute2DCoords(m2)
        molecule_2d = AllChem.EmbedMolecule(m2)

        return molecule_2d
    
    def convert3D(self, molecule_2d):
        """Convert 2D geometry into 3D geometry using RDKit"""

        molecule_3d = AllChem.MMFFOptimizeMolecule(molecule_2d)

        return molecule_3d
    
    def save(self, molecule_3d):
        """Write 3D molecule to file (xyz and mol/mol2)"""

        Chem.MoltoMolBlock(m2)
        file=open('path/file.mol', 'w+') 

    def total_partial_charge(self):
        return np.array([a.partialcharge for a in self.atoms]).sum()

    def natoms(self):
        return len(self.atoms)

    def pop_atom(path, output, atom='Na'):
        to_remove = []
        to_save = []
        with open(path, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if i == 2:
                    info = [int(x) for x in line.split()]
                elif atom.upper() in line:
                    to_remove.append(i)
                    to_save.append(line)

        change = len(to_remove)
        with open(output, 'w') as f:
            for i, line in enumerate(lines):
                if i == 2:
                    info[0] -= change
                    f.write(' %s %s %s %s %s\n' % tuple(info))
                elif i in to_remove:
                    pass
                else:
                    f.write(line)

        return to_remove, to_save

    def push_atom(path, output, idx, content):
        with open(path, 'r') as f:
            lines = f.readlines()

        info = [int(x) for x in lines[2].split()]

        for i, line in zip(idx, content):
            if 'NA' in line:
                parts = line.split()
                parts[1] = 'NA'
                parts[5] = 'Na+'
                parts[6] = '2'
                parts[7] = 'Na+'
                parts[8] = '1.0000'
                line = ' '.join(parts) + '\n'
            lines.insert(i + 1, line)

        change = len(idx)
        with open(output, 'w') as f:
            for i, line in enumerate(lines):
                if i == 2:
                    info[0] += change
                    f.write(' %s %s %s %s %s\n' % tuple(info))
                else:
                    f.write(line)


