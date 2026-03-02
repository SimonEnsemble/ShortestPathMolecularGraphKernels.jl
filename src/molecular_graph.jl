struct MolGraph
    smiles::String
    incl_hydrogen::Bool
    g::SimpleGraph
    atoms::Vector{AtomType}
    bonds::Vector{BondType}
    spaths::Vector{ShortestPath}
end

function MolGraph(smiles::String; incl_hydrogen::Bool=false, verbose::Bool=false)
    if verbose
        println("reading in ", smiles)
        println("\tignore H atoms? ", ! incl_hydrogen)
    end

    # use MolecularGraph.jl to read SMILES
    mol = smilestomol(smiles)
    if incl_hydrogen
        add_hydrogens!(mol)
    end

    # construct simple graph
    g = deepcopy(mol.graph)
    if ! is_connected(g)
        @warn "molecular graph for $smiles is not connected!"
    end
    
    # node labels: atom type
    atoms = [AtomType(elements[mol.vprops[v].symbol].number) for v = 1:nv(g)]
    
    # edge labels: bond order
    aromatic_bond = is_edge_aromatic(mol)
    bonds = [
        aromatic_bond[i] ? BondType(0) : BondType(mol.eprops[ed].order)
            for (i, ed) in enumerate(edges(g))
    ]

    # shortest path lengths between all pairs of nodes
    #   (lazy computation; don't do this unless called.)
    spaths = ShortestPath[]
    
    return MolGraph(smiles, incl_hydrogen, g, atoms, bonds, spaths)
end

get_shortest_paths(mg::MolGraph, u::Int, v::Int) = get_shortest_paths(mg.spaths, u, v)
