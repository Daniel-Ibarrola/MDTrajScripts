

def extract_chain(traj, chains):
    """ Extract chains from a trajectory

        Parameters
        ----------
        traj: mdtraj.trajectory

        chains: list[int]

        Returns
        -------
        new_traj: mdtraj.Trajectory

    """
    topology = traj.topology
    return traj.atom_slice(
        [atom.index for atom in topology.atoms if (atom.residue.chain.index in chains)]
    )
