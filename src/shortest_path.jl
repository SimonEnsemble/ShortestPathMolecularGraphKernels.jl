"""
    ShortestPath(
        ρ::Vector{Int}, seq::Vector{UInt8}, ℓ::Int, w::Float64
    ) -> ShortestPath

a data structure representing a shortest path between two nodes in a molecular graph.


### fields
* `ρ::Vector{Int}`: the path represented as an ordered sequence of node indices
* `seq::Vector{UInt8}`: the atom-bond label sequence traversed along the path
* `ℓ::Int`: the length of the path (= # of hops)
* `w::Float64`: path weight, which ≂̸ 1.0 only when multiple degenerate 
shortest paths exist between the same pair of nodes, used to normalize contributions.

_note_: under the hood, we store the shortest paths in a canonical order to avoid 
compute-intensive reversal for matching.
"""
struct ShortestPath
	ρ::Vector{Int}
    seq::Vector{UInt8}
	ℓ::Int
	w::Float64
end

function Base.show(io::IO, sp::ShortestPath)
    println(io, "length: ", sp.ℓ)
    println(io, "weight: ", sp.w)
    println(io, "ρ = ", sp.ρ)
    println(io, "atom sequence: ", 
            [AtomType(sp.seq[2*i-1]) for i = 1:sp.ℓ+1]
    )
    println(io, "bond sequence: ", 
            [BondType(sp.seq[2*i]) for i = 1:sp.ℓ]
    )
end

"""
    get_shortest_paths(spaths::Vector{ShortestPath}, u::Int, v::Int) -> Vector{ShortestPath}
    get_shortest_paths(mg::MolGraph, u::Int, v::Int) -> Vector{ShortestPath}
    get_shortest_paths(mg::MolGraph, u::Int) -> Vector{ShortestPath}

filter the list of shortest paths to get those "u ↔ v" or "u ↔ *"
"""
function get_shortest_paths(spaths::Vector{ShortestPath}, u::Int, v::Int)
    return filter(spaths) do spath
        ρ_first = spath.ρ[1]
        ρ_last  = spath.ρ[end]
        return (u == ρ_first && v == ρ_last) || (v == ρ_first && u == ρ_last)
    end
end

function get_shortest_paths(spaths::Vector{ShortestPath}, u::Int)
    return filter(spaths) do spath
        return u == spath.ρ[1] || u == spath.ρ[end]
    end
end
