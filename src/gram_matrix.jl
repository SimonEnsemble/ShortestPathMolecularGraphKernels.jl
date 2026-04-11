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

"""
    center_Gram_matrix(K::Matrix{Float64}) -> Matrix{Float64}

center a Gram (kernel) matrix using double centering.

double centering of a Gram matrix removes the mean along each row and each column, 
then adds back the grand mean to avoid over-subtraction. this is equivalent to 
re-computing the original feature vectors, centering them, *then* taking the dot product 
of those centered feature vectors. do this before e.g. k-PCA!

# example
```julia
Xs = rand(2, 10) # 10 2D feature vectors
K = Xs' * Xs
K_centered = center_Gram_matrix(K)
```
"""
function center_Gram_matrix(K::Matrix{Float64})
    @assert size(K, 2) == size(K, 1) "Gram matrix must be square!"

    row_mean = mean(K, dims=2)
    col_mean = mean(K, dims=1)
    total_mean = mean(K)

    return K .- row_mean .- col_mean .+ total_mean
end
