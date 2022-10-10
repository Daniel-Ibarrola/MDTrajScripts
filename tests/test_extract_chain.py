from ..src.extract_chain import extract_chain
import mdtraj as mdt


def test_extract_single_chain():
    traj = mdt.load("../test_cases/eralpha/1qku/1qku.pdb")

    new_traj = extract_chain(traj, [0])
    assert new_traj.n_chains == 1


def test_extract_multiple_chains():
    traj = mdt.load("../test_cases/eralpha/1qku/1qku.pdb")

    new_traj = extract_chain(traj, [0, 1, 2])
    assert new_traj.n_chains == 3
