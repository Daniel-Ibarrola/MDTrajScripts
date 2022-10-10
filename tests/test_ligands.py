import pytest
import numpy as np
from ..src.ligands import (find_ligands_in_traj, get_ligand_atom_indices, ligand_centroid,
maximum_distance, ligand_maximum_extent)
import mdtraj as mdt


@pytest.fixture
def trajectory():
    return mdt.load("../test_cases/eralpha/1qku/1qku.pdb")


def test_find_ligands(trajectory):
    ligands = find_ligands_in_traj(trajectory)
    assert len(ligands) == 3
    assert ligands == [
        "EST:D",
        "EST:E",
        "EST:F",
    ]


def test_get_ligand_atom_indices(trajectory):
    traj = trajectory
    ligand_indices = get_ligand_atom_indices(traj, "EST:D")
    expected_indices = list(range(5940, 5960))
    assert ligand_indices == expected_indices


def test_ligand_centroid(trajectory):
    traj = trajectory
    ligand_id = "EST:D"
    centroid = ligand_centroid(traj, ligand_id)

    coordinates = np.array(
        [[[10.4106, 1.7203, 2.4775],
          [10.2995, 1.7834, 2.537],
          [10.1695, 1.7355, 2.512],
          [10.0598, 1.799, 2.5704],
          [10.1506, 1.624, 2.4274],
          [10.2621, 1.5588, 2.366],
          [10.2371, 1.4379, 2.2735],
          [10.3644, 1.3753, 2.2086],
          [10.4898, 1.3873, 2.2953],
          [10.5178, 1.5388, 2.3261],
          [10.3957, 1.6078, 2.3918],
          [10.6462, 1.5459, 2.4125],
          [10.7711, 1.4803, 2.3508],
          [10.7463, 1.3343, 2.3124],
          [10.617, 1.327, 2.2242],
          [10.6228, 1.1821, 2.1792],
          [10.7701, 1.1713, 2.1263],
          [10.8494, 1.2719, 2.2135],
          [10.961, 1.2027, 2.2746],
          [10.7379, 1.2449, 2.4419]]],
        dtype=np.float32)
    expected_centroid = np.mean(coordinates, axis=1)[0]
    assert centroid.shape == (3, )
    assert np.allclose(centroid, expected_centroid)


def test_maximum_distance():
    centroid = np.array([0., 0., 0.])
    coordinates = np.array([
        [1., 1., 1.],
        [2., 2., 2.],
        [3., 3., 3.],
        [4., 4., 4.],
    ])
    assert maximum_distance(centroid, coordinates) == np.sqrt(48)


def test_maximum_ligand_extent(trajectory):
    assert np.allclose(ligand_maximum_extent(trajectory, "EST:D"), 0.59849)
