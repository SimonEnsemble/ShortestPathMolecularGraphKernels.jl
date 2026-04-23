function get_edge_id(g::SimpleGraph, u::Int, v::Int)
	for (i, ed) in enumerate(edges(g))
		if (ed.src == u && ed.dst == v) || (ed.src == v && ed.dst == u)
			return i
		end
	end
    error("no such edge ($u, $v) in the graph!")
end

function nb_spaths(mg::MolGraph)
	p = 0 # number of shortest paths
	for cc in connected_components(mg.g)
		# paths between nodes in this connected component
		n = length(cc)
		p += n * (n - 1) ÷ 2 
	end
	return p + nv(mg.g) # self-loops
end

"""
    find_shortest_paths!(mg::MolGraph; ℓ_max::Int=typemax(Int)) -> nothing

search for all shortest paths between all pairs of nodes in a molecular graph `mg`,
up to and including of length `ℓ_max`.
(uses Dijkstra's algorithm.)
attach the list of shortest paths and their atom-bond label sequences to the molecular graph.

# example
```julia
caf = MolGraph("CN1C=NC2=C1C(=O)N(C)C(=O)N2C") # caffeine
caf.spaths # empty array
find_shortest_paths!(caf)
caf.spaths # now contains the list of shortest paths
```
"""
function find_shortest_paths!(mg::MolGraph; ℓ_max::Int=typemax(Int))
	# loop over (src, dst) pairs
	for src = 1:nv(mg.g)
		# find all (src, ..., dst) paths for ∀ dst
		dsp = dijkstra_shortest_paths(mg.g, src, allpaths=true)

		for dst = src:nv(mg.g)
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
                    seq[2 * i - 1] = UInt8(mg.atoms[u])
                    
                    if i != ℓ+1
                        v = path[i+1]
                        i_edge = get_edge_id(mg.g, u, v)

                        # note bond label
                        seq[2 * i] = UInt8(mg.bonds[i_edge])
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

                push!(mg.spaths, sp)
			end			
		end
	end

    @assert sum(sp.w for sp in mg.spaths) ≈ nb_spaths(mg) * 1.0

    return nothing
end
