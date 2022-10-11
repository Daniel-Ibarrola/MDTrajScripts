import numpy as np

# Distance from the ligand centroid to the binding site
# value from plip
BS_DIST = 0.85  # in nanometers


def get_binding_site_atoms_indices(lig_center, lig_max_extent, coords):
    """ Get the indices of all the atoms that belong to the binding site.

        Parameters
        ----------
        lig_center : np.ndarray of shape (3,)
        lig_max_extent : float
        coords: np.ndarray of shape(n_atoms, 3)

        Returns
        -------
        bs_indices: np.ndarray
    """
    bs_cutoff = lig_max_extent + BS_DIST
    distance = np.sqrt(np.sum(np.power(coords - lig_center, 2), axis=1))
    return np.where(np.logical_and(distance > lig_max_extent, distance <= bs_cutoff))[0]
