import mdtraj as mdt


non_ligand_residues = [
    "ALA", "ARG", "ASN", "ASP", "ASX",
    "CYS", "GLU", "GLN", "GLX", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HOH", "DOD",
]


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

    traj = mdt.load(topology_file)
    topology = traj.topology
    n_residues = topology.n_residues

    # Find all possible ligand residues
    ligands = []
    for ii in range(n_residues):
        residue = topology.residue(ii)
        if residue.name not in non_ligand_residues \
                and residue.n_atoms > 5 \
                and residue.name not in ligands:
            ligands.append(residue.name)

    return ligands
