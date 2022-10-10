from ..src.ligands import find_ligands_in_traj, get_ligand_atom_indices
import os
import mdtraj as mdt


def test_can_load_file():
    assert os.path.isfile("../test_cases/eralpha/1qku/1qku.pdb")


def test_find_ligands():
    ligands = find_ligands_in_traj("../test_cases/eralpha/1qku/1qku.pdb")
    assert len(ligands) == 3
    assert ligands == [
        "EST:D",
        "EST:E",
        "EST:F",
    ]


def test_get_ligand_atom_indices():
    traj = mdt.load("../test_cases/eralpha/1qku/1qku.pdb")
    ligand_indices = get_ligand_atom_indices(traj, "EST:D")
    expected_indices = list(range(5940, 5960))
    assert ligand_indices == expected_indices
