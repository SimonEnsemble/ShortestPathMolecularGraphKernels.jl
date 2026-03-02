struct ShortestPath
    # path (sequence of nodes) in the graph
	ρ::Vector{Int}
    # atom-edge label sequence
    seq::Vector{UInt8}
    # length
	ℓ::Int
    # weight (≠1 when multiple shortest paths exist between a pair of nodes)
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

function get_shortest_paths(spaths::Vector{ShortestPath}, u::Int, v::Int)
    return filter(spaths) do spath
        ρ_first = spath.ρ[1]
        ρ_last  = spath.ρ[end]
        return (u == ρ_first && v == ρ_last) || (v == ρ_first && u == ρ_last)
    end
end
