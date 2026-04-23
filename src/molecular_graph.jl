struct MolGraph
    smiles::String
    incl_hydrogen::Bool
    g::SimpleGraph
    atoms::Vector{AtomType}
    bonds::Vector{BondType}
    spaths::Vector{ShortestPath}
end

"""
    MolGraph(
        smiles::String; incl_hydrogen::Bool=false, verbose::Bool=false
    ) -> MolGraph

interpret a SMILES representation of a molecular structure as a molecular graph.

relies on `MolecularGraph.jl`.

# fields
- `smiles::String`: the SMILES of the molecule represented
- `incl_hydrogen::Bool`: include H atoms?
- `g::SimpleGraph`: the graph (nodes: atoms; edges: bonds)
- `atoms::Vector{AtomType}`: list of labels (atom types) on the nodes
- `bonds::Vector{BondType}`: list of labels (bond types) on the edges
- `spaths::Vector{ShortestPath}`: list of shortest paths

# example
```julia
mg = MolGraph("CN1C=NC2=C1C(=O)N(C)C(=O)N2C") # caffeine
```
"""
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
get_shortest_paths(mg::MolGraph, u::Int) = get_shortest_paths(mg.spaths, u)
