__precompile__(false)
@warn "precompile set to false"

module ShortestPathMolecularGraphKernels

import Base: ==, hash
using Base.Threads

using Graphs, CairoMakie, GraphMakie, Printf, Colors
using PeriodicTable: elements
using MolecularGraph: smilestomol, SMILESMolGraph, is_edge_aromatic, 
    coords2d, RASMOL_ATOM_COLOR, coordgen!, add_hydrogens!
set_theme!(
    Axis = (
        xautolimitmargin = (0.1, 0.1),
        yautolimitmargin = (0.1, 0.1),
    )
)

# atom and bond encodings (export all for use in testing)
include("atom_bond_encodings.jl")

export BondType
for instance in instances(BondType)
    @eval export $(Symbol(instance))
end

export AtomType
for instance in instances(AtomType)
    @eval export $(Symbol(instance))
end

include("shortest_path.jl")
include("molecular_graph.jl")
include("viz.jl")
include("find_shortest_paths.jl")
include("shortest_path_graph_kernel.jl")
include("featurization.jl")
include("gram_matrix.jl")

export ShortestPath, get_shortest_paths,                          # shortest_path.jl 
    MolGraph,                                                     # molecular_graph.jl
    viz,                                                          # viz.jl
    find_shortest_paths!,                                         # find_shortest_paths.jl
    OrdinarySPF, LabelSeqSPF, build_feature_vector, FeatureSpace, # featurization.jl
    active_features,
    shortest_path_graph_kernel, shortest_path_kernel, indistinguishable_pairs, # shorest_path_kernel.jl
    compute_Gram_matrix # gram_matrix.jl
end
