import mdtraj as mdt


chain_names = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]


def is_ligand_atom(atom):
    if not atom.residue.is_water and not atom.residue.is_protein\
            and atom.residue.n_atoms > 4:
        return True
    return False


def find_ligands_in_traj(topology_file):
    """ Returns a list of ligand ids in a trajectory.

        Parameters
        ----------
        topology_file : str
            A topology file such as a pdb.

        Returns
        -------
        ligands : list[str]
            A list of the ligands ids in the topology file
    """
    topology: mdt.Topology
    traj: mdt.Trajectory

    traj = mdt.load(topology_file)
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
            A topology file such as a pdb.

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
