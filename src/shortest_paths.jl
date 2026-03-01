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

function get_edge_id(g::SimpleGraph, u::Int, v::Int)
	for (i, ed) in enumerate(edges(g))
		if (ed.src == u && ed.dst == v) || (ed.src == v && ed.dst == u)
			return i
		end
	end
    error("no such edge ($u, $v) in the graph!")
end

function get_shortest_paths(spaths::Vector{ShortestPath}, u::Int, v::Int)
    return filter(spaths) do spath
        ρ_first = spath.ρ[1]
        ρ_last  = spath.ρ[end]
        return (u == ρ_first && v == ρ_last) || (v == ρ_first && u == ρ_last)
    end
end

function find_all_shortest_paths(
    g::SimpleGraph, atoms::Vector{AtomType}, bonds::Vector{BondType}; ℓ_max::Int=typemax(Int)
)
    @assert is_connected(g)

	sps = ShortestPath[]
	# loop over (src, dst) pairs
	for src = 1:nv(g)
		# find all (src, ..., dst) paths for ∀ dst
		dsp = dijkstra_shortest_paths(g, src, allpaths=true)

		for dst = src:nv(g)
			#=
			construct all (src, ..., dst) paths
			=#
			all_paths = Vector{Vector{Int}}()
			
			stack = [(dst, Int[])]
		    while ! isempty(stack)
		        current_node, current_path = pop!(stack)
		        
		        # append current node (building backward: dst → ... → src)
		        new_path = push!(copy(current_path), current_node)
		        
		        if current_node == src # finished!
		            push!(all_paths, reverse(new_path))
		            continue
		        end
		        
		        for pred in dsp.predecessors[current_node] # not finished.
		            push!(stack, (pred, new_path))
		        end
            end

			@assert length(all_paths) == Int(dsp.pathcounts[dst])

			#=
			store paths we found
			=#
			for path in all_paths
                # down-weigh multiple shortest paths between the same node pair
				w = 1.0 / dsp.pathcounts[dst]

                # length (# hops = # edges)
                ℓ = length(path) - 1
                if ℓ > ℓ_max
                    continue # do not store
                end

                #=
                build atom-bond label sequence
                =#
                seq = Vector{UInt8}(undef, 2*ℓ+1)
                for i = 1:ℓ+1 # loop over nodes
                    # looking at hop (u, v)
                    u = path[i]

                    # note atom label
                    seq[2 * i - 1] = UInt8(atoms[u])
                    
                    if i != ℓ+1
                        v = path[i+1]
                        i_edge = get_edge_id(g, u, v)

                        # note bond label
                        seq[2 * i] = UInt8(bonds[i_edge])
                    end   
                end

                #=
                note path, but first canonicalize for easy comparison
                =#
                reverse_seq = reverse(seq)
                if isless(seq, reverse_seq) # atoms in a canonical order
                    sp = ShortestPath(reverse(path), reverse_seq, ℓ, w) # good as is
                else
                    sp = ShortestPath(path, seq, ℓ, w) # good as is
                end

                push!(sps, sp)
			end			
		end
	end

	# only if connected...
    @assert sum(sp.w for sp in sps) ≈ nv(g) * (nv(g) - 1) / 2 + nv(g)
	
	return sps
end
