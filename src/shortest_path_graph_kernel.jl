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

"""
    shortest_path_graph_kernel(
        mg₁::MolGraph, mg₂::MolGraph; exact_seq_matching::Bool=false
    ) -> Float64

compute the shortest path molecular graph kernel between two input
molecular graphs.

the kernel value is calculated by summing the similarity of all pairs of shortest
paths between the two graphs. 

if `exact_seq_matching::Bool`, we require exact matching of the atom-bond label
sequence to note two common paths; otherwise, we note two common paths if the
length and unordered pair of atom labels in the src and dst nodes are identical.

returns a numerical value representing the similarity between the two graphs.
"""
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
