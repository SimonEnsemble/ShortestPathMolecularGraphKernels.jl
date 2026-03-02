"""
    compute_Gram_matrix(mgs::Vector{MolGraph}, exact_seq_matching::Bool) -> Matrix{Float64}

compute the symmetric Gram (kernel) matrix for a list of molecular graphs.

Each element `K[i, j]` in the resulting matrix `K` represents the similarity between 
molecular graph `i` and `j` according to the shortest path graph kernel.

### implementation details
* **symmetry:** to reduce computation time by 50%, only the upper triangle is 
  explicitly calculated, and values are mirrored to the lower triangle.
* **parallelization:** this function utilizes `Base.Threads.@threads` to distribute 
  the outer loop iterations across available CPU cores.
"""
function compute_Gram_matrix(
	mgs::Vector{MolGraph}, exact_seq_matching::Bool
)
	n = length(mgs)
	K = zeros(length(mgs), length(mgs))
	@threads for i = 1:length(mgs)
		for j = i:length(mgs)
			K[i, j] = K[j, i] = shortest_path_graph_kernel(
				mgs[i], mgs[j], exact_seq_matching=exact_seq_matching
			)
		end
	end
	return K
end
