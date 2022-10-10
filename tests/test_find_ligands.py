from ..src import find_ligands
import os


def test_can_load_file():
    assert os.path.isfile("../test_cases/eralpha/1qku/1qku.pdb")


def test_find_ligands():
    ligands = find_ligands.find_ligands_in_traj("../test_cases/eralpha/1qku/1qku.pdb")
    assert len(ligands) == 3

