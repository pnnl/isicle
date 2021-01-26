import pytest
import isicle
import os
import pandas as pd

def localfile(path):
    "Returns path relative to this file."
    return os.path.join(os.path.dirname(__file__), path)

@pytest.fixture()
def moleculeprep(path):
    return isicle.molecule_prep.MolecularString()

class TestMolecularString:

    def test_init(self, moleculeprep, path):
        assert isinstance(moleculeprep, isicle.molecule_prep.MolecularString)

    @pytest.mark.parametrize('path,expected',
                             [('resources/sample_smiles_plain.smi', "CCOP(C)(=O)SCCN(C(C)C)C(C)C"),
                              ('resources/sample_no_smiles.smi', None)])
    def test_load(self, moleculeprep, path, expected):
        # initialize
        contents = moleculeprep.load(localfile(path))

        # test attribute
        assert moleculeprep.contents == expected

        # test return
        assert contents == expected

    @pytest.mark.parametrize('path,expected',
                             [('resources/sample_smiles_desalt.smi', "CCOP(C)(=O)SCCN(C(C)C)C(C)C"),
                              ('resources/sample_smiles_plain.smi', "CCOP(C)(=O)SCCN(C(C)C)C(C)C")])
    def test_desalt(self, moleculeprep, path, expected):
        # initialize
        moleculeprep.load(localfile(path))
        result = moleculeprep.desalt()

        # test attribute
        assert moleculeprep.result == expected

        # test return
        assert result == expected

    @pytest.mark.parametrize('path,expected',
                             [('resources/sample_smiles_neutral.smi', "CCOP(C)(=O)SCCN(C(C)C)C(C)C"),
                              ('resources/sample_smiles_plain.smi', "CCOP(C)(=O)SCCN(C(C)C)C(C)C")])
    def test_neutralize(self, moleculeprep, path, expected):
        # initialize
        moleculeprep.load(localfile(path))
        result = moleculeprep.neutralize()

        # test attribute
        assert moleculeprep.result == expected

        # test return
        assert result == expected

    @pytest.mark.parametrize('path,expected',
                             ('resources/sample_smiles_plain.smi', "CCOP(C)(=O)SCCN(C(C)C)C(C)C"))
    def test_tautomerize(self, moleculeprep, path, expected):
        # initialize
        moleculeprep.load(localfile(path))
        result = moleculeprep.tautomerize()

        # test attribute
        assert moleculeprep.result == expected

        # test return
        assert result == expected

"""
    @pytest.mark.parametrize('path,expected',
                             [('resources/sample_mol.mol',1)
                              ('resources/sample_no_mol.mol',0)])
    def test_to_smiles(self,moleculeprep,path,expected):
        #initialize
        moleculeprep.load(localfile(path))
        result = moleculeprep.to_smiles()

        # test attribute
        assert len(moleculeprep.result) == expected

        # test return
        assert len(result) == expected


    @pytest.mark.parametrize('path,output,expected',
                             ('resources/sample_smiles_plain.smi', 'resources/sample_smi.out',"CCOP(C)(=O)SCCN(C(C)C)C(C)C"))
    def test_to_file(self,moleculeprep,path,expected):
        #initialize
        moleculeprep.load(localfile(path))
        moleculeprep.to_file(type=string,output)
        with open(output, 'r') as f:
            result = f.readlines()

        #test attribute
        assert moleculeprep.to_file

        #test return
        assert result == expected
"""
