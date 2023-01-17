from ..src.to_rdkit import to_rdkit
import mdtraj
import pytest


@pytest.fixture
def estradiol():
    return mdtraj.load("data/estradiol.pdb")


def test_traj_to_rdkit(estradiol):
    mol = to_rdkit(estradiol)
    assert mol.GetNumAtoms() == 20

