function compute_Gram_matrix(
	mgs::Vector{MolGraph}, exact_seq_matching::Bool
)
	n = length(mgs)
	K = zeros(length(mgs), length(mgs))
	for i = 1:length(mgs) # TODO parallelize
		for j = i:length(mgs)
			K[i, j] = K[j, i] = shortest_path_graph_kernel(
				mgs[i], mgs[j], exact_seq_matching=exact_seq_matching
			)
		end
	end
	return K
end
