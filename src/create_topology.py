import mdtraj as mdt
import mdtraj.core.element as mdt_element


def create_topology():
    """ Create a Topology object from scratch.
    """
    topology = mdt.Topology()

    # Adding a new chain is pretty straightforward
    chain = topology.add_chain()
    # Adding a residue
    residue = topology.add_residue(name="ASP", chain=chain, resSeq=1)
    # Adding an atom
    element = mdt_element.get_by_symbol("C")
    topology.add_atom(name="CA",
                      element=element,
                      residue=residue)

    element = mdt_element.get_by_symbol("N")
    topology.add_atom(name="NZ",
                      element=element,
                      residue=residue)

    # Adding a bond
    topology.add_bond(atom1=topology.atom(0),
                      atom2=topology.atom(1)
                      )
    return topology
