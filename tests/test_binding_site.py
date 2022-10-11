from ..src.binding_site import get_binding_site_atoms_indices
import numpy as np


def test_get_binding_site_atoms_indices():
    ligand_maximum_extent = 0.15
    ligand_centroid = np.array([0.0, 0.0, 0.0])
    coordinates = np.array([
        [0.01, 0.01, 0.01],
        [0.02, 0.02, 0.02],
        [0.3, 0.3, 0.3],
        [0.4, 0.4, 0.4],
        [0.5, 0.5, 0.5],
        [1., 1., 1.],
        [2., 2., 2.],
    ])
    bsite_indices = get_binding_site_atoms_indices(
        ligand_centroid, ligand_maximum_extent, coordinates)
    assert isinstance(bsite_indices, np.ndarray)

    expected_indices = np.array([2, 3, 4])
    assert np.all(bsite_indices == expected_indices)
