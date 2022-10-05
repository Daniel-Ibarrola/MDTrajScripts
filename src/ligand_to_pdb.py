import mdtraj as mdt
import os
from rdkit import Chem


non_ligand_residues = [
    "ALA", "ARG", "ASN", "ASP", "ASX",
    "CYS", "GLU", "GLN", "GLX", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HOH", "DOD",
]


def residue_to_rdkit(residue_index, topology):
    """ Transform a residue into a rdkit molecule. """
    topology: mdt.Topology

    molecule = Chem.RWMol()
    residue_atoms = []
    atom_index_mapper = {}

    for atom in topology.residue(residue_index).atoms:
        index_topology = atom.index
        element = atom.element.symbol
        if element == "D":
            element = "H"

        index_rd_mol = molecule.AddAtom(Chem.Atom(element))
        residue_atoms.append(index_topology)
        atom_index_mapper[index_topology] = index_rd_mol

    bonds = []
    for atom_1, atom_2 in topology.bonds:
        if atom_1.index in residue_atoms and atom_2.index in residue_atoms:
            index_1 = atom_index_mapper[atom_1.index]
            index_2 = atom_index_mapper[atom_2.index]
            if (index_1, index_2) not in bonds:
                molecule.AddBond(index_1, index_2)
                bonds.append((index_1, index_2))

    return molecule


def ligands_to_pdb(pdb_file, save_to):
    """ Given a pdb file, extracts the ligand and saves it to pdb.

        Parameters
        ----------
        pdb_file : str
            Path to a pdb file.

        save_to : str
            The path where the extracted ligands will be saved.
    """
    topology: mdt.Topology

    traj = mdt.load(pdb_file)
    topology = traj.topology
    ligands_residues = []  # index of the residues

    for ii in range(topology.n_residues):
        residue = topology.residue(ii)
        if residue.name not in non_ligand_residues \
                and residue.n_atoms > 5 \
                and ii not in ligands_residues:
            ligands_residues.append(ii)

    for residue in ligands_residues:
        ligand_name = topology.residue(residue).name
        try:
            molecule = residue_to_rdkit(residue, topology)
        except ValueError:
            print(f"Error in {pdb_file} with ligand {ligand_name}")
            raise ValueError
        pdb_string = Chem.MolToPDBBlock(molecule)
        with open(os.path.join(save_to, ligand_name + ".pdb"), "w") as fp:
            fp.write(pdb_string)
