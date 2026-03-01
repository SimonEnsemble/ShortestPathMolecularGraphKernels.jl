#=
define shortest path features
=#
abstract type ShortestPathFeature end

struct OrdinarySPF <: ShortestPathFeature
    lᵤ::UInt8
    lᵥ::UInt8
    ℓ::Int
end
OrdinarySPF(ℓᵤ::AtomType, ℓᵥ::AtomType, ℓ::Int) = OrdinarySPF(
    UInt8(ℓᵤ), UInt8(ℓᵥ), ℓ
)

struct LabelSeqSPF <: ShortestPathFeature
    seq::Vector{UInt8}
end

function ==(f1::LabelSeqSPF, f2::LabelSeqSPF)
    return f1.seq == f2.seq
end

function hash(f::LabelSeqSPF, h::UInt)
    return hash((:LabelSeqSPF, f.seq), h)
end

#=
distill shortest path into the feature it activates
=#
OrdinarySPF(sp::ShortestPath) = OrdinarySPF(
    minmax(sp.seq[1], sp.seq[end])..., sp.ℓ
)

LabelSeqSPF(sp::ShortestPath) = LabelSeqSPF(sp.seq)

#=
obtain list of activated features in a list of molecular graphs
=#
function active_features(mgs::Vector{MolGraph}, ::Type{T}) where {
	T <: ShortestPathFeature}
	found_features = Set{T}()
	for mg in mgs 
		for sp in mg.spaths
			push!(found_features, T(sp))
		end
	end
	return collect(found_features)
end

#=
construct feature space based on activate features in a set of molecular graphs
=#
struct FeatureSpace{T <: ShortestPathFeature}
    features::Vector{T}
    feature_to_id::Dict{T, Int}
    dims::Int
end

function FeatureSpace{T}(mgs::Vector{MolGraph}) where {
    T <: ShortestPathFeature}
    features = active_features(mgs, T)
    
    feature_to_id = Dict{T, Int}(f => i for (i, f) in enumerate(features))

    return FeatureSpace{T}(features, feature_to_id, length(features))
end

#=
build feature vector
=#
function build_feature_vector(mg::MolGraph, fs::FeatureSpace{T}) where {
	T <: ShortestPathFeature
}
    ϕ = zeros(fs.dims)
    for sp in mg.spaths
        feature = T(sp)
        id = fs.feature_to_id[feature]
        ϕ[id] += sp.w
    end
    return ϕ
end
