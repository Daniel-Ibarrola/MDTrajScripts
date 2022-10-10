import mdtraj as mdt
import numpy as np


chain_names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]


def is_ligand_atom(atom):
    """ Check if an atom belongs to a ligand. """
    if not atom.residue.is_water and not atom.residue.is_protein\
            and atom.residue.n_atoms > 4:
        return True
    return False


def find_ligands_in_traj(traj):
    """ Returns a list of ligand ids in a trajectory.

        Parameters
        ----------
        traj : mdtraj.Trajectory

        Returns
        -------
        ligands : list[str]
            A list of the ligands ids in the topology file
    """
    topology = traj.topology

    ligands = []
    for atom in topology.atoms:
        if is_ligand_atom(atom):
            chain = atom.residue.chain.index
            ligand_id = atom.residue.name + ":" + chain_names[chain]
            if ligand_id not in ligands:
                ligands.append(ligand_id)

    return ligands


def get_ligand_atom_indices(traj, ligand_id):
    """ Returns the indices of the atoms of the ligand with the given id
        in the trajectory

        Parameters
        ----------
        traj : mdtraj.Trajectory

        ligand_id : str

        Returns
        -------
        indices : list[int]

    """
    topology = traj.topology
    indices = []

    ligand, chain = ligand_id.split(":")
    chain_index = chain_names.index(chain)
    for atom in topology.atoms:
        if atom.residue.name == ligand and atom.residue.chain.index == chain_index:
            indices.append(atom.index)

    return indices


def ligand_centroid(traj, ligand_id):
    """ Get the centroid of the ligand with the given id."""
    indices = get_ligand_atom_indices(traj, ligand_id)
    coords = traj.xyz[:, indices, :]
    return np.mean(coords, axis=1)[0]


def maximum_distance(centroid, coordinates):
    """ Get the maximum distance from the centroid to the given
        coordinates
    """
    distance = np.sqrt(np.sum(np.power(coordinates - centroid, 2), axis=1))
    return np.amax(distance)


def ligand_maximum_extent(traj, ligand_id):
    """ Computes the maximum extent of the ligand. This is the maximum distance
        from the centroid to any of its atoms.
    """
    centroid = ligand_centroid(traj, ligand_id)
    indices = get_ligand_atom_indices(traj, ligand_id)
    coords = traj.xyz[:, indices, :][0]
    return maximum_distance(centroid, coords)
