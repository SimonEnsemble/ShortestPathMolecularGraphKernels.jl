function same_src_dst_label_pair(p₁::ShortestPath, p₂::ShortestPath)
    return minmax(p₁.seq[1], p₁.seq[end]) == minmax(p₂.seq[1], p₂.seq[end])
end

function shortest_path_kernel(
    p₁::ShortestPath, p₂::ShortestPath, exact_seq_matching::Bool
)
    k = 0.0
    # same length?
    if p₁.ℓ == p₂.ℓ
        if exact_seq_matching
            # canonicalized, so don't need the reverse direction
            if p₁.seq == p₂.seq
                k += p₁.w * p₂.w
            end
        else
            if same_src_dst_label_pair(p₁, p₂)
                k += p₁.w * p₂.w
            end
        end
    end
    # no similarity
    return k
end

function shortest_path_graph_kernel(
    mg₁::MolGraph, mg₂::MolGraph; exact_seq_matching::Bool=false
)
	k = 0 # will be kernel value
    # loop over shortest paths in molecular graph #1
    for sp₁ in mg₁.spaths
        # loop over shortest paths in molecular graph #2
        for sp₂ in mg₂.spaths
            # increment molecular graph similarity with shortest-path similarity
            k += shortest_path_kernel(sp₁, sp₂, exact_seq_matching)
        end
    end
	
	return k
end

function indistinguishable_pairs(K::Matrix{Float64})
	n = size(K)[1] # number of molecules
	psps = Vector{Tuple{Int, Int}}() # list of perfectly similar pairs
	# loop thru pairs
	for i = 1:n
		for j = i+1:n
			if all(isapprox.(K[i, :], K[j, :]))
				push!(psps, (i, j))
			end
		end
	end
	return psps
end
