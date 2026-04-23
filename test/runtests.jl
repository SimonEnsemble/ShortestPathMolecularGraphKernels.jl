### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ cecf3058-bb8f-11f0-97f3-bda46249b7c9
begin
	# if running as a Pluto notebook...
	if isdefined(Main, :PlutoRunner)
		import Pkg; Pkg.activate("../")
		using Revise
	end
	
	using ShortestPathMolecularGraphKernels
	using PlutoUI # sorry no way around for runtests.jl as script
	using Graphs, Printf, LinearAlgebra, GraphMakie, BenchmarkTools, Colors, Test

	# if u want a cool theme
	# using MakieThemes, CairoMakie
	# set_theme!(ggthemr(:pale))

	TableOfContents()
end

# ╔═╡ f9461018-9522-4fe6-a53f-665b3839b6e4
md"# `ShortestPathMolecularGraphKernels.jl`: demo and tests

## 🐧 constructing the representation of a molecule

1. interpret a SMILES string as a molecular graph (nodes: atoms; edges: labels; categorical node attribute: atom type; categorical edge attribute: bond type)
2. search for all shortest paths between all pairs of nodes and store them---including the atom-bond label sequence along them
3. visualize the molecular graph and explore its shortest paths
"

# ╔═╡ 50d85bc2-2e1a-424a-84f7-5a22c88759ed
begin
	# interpret the SMILES string as a molecular graph
	mg = MolGraph(
		"CN1C=NC2=C1C(=O)N(C)C(=O)N2C", # SMILES
		incl_hydrogen=false             # "include H atoms?
	)
	
	# search for and store all shortest paths---
	#   including the atom and bond label sequence
	find_shortest_paths!(
		mg,                 # the molecular graph
		ℓ_max=typemax(Int)  # the max path length
	)

	# show the attributes of the data structure
	mg
end

# ╔═╡ f774c2df-1859-443e-9fee-9d490a18c83b
md"id of shortest path $(@bind id_spath PlutoUI.Slider(1:length(mg.spaths)))"

# ╔═╡ e978a83d-f0b3-4dd0-98e6-250124698ec4
id_spath

# ╔═╡ ffa9d1c9-9534-46fd-aed2-ed549ac46384
mg.spaths[id_spath] # wut path are we lookin' at?

# ╔═╡ 1f17de1a-5375-4c22-97cd-7dc04645b9d1
viz(mg, id_spath=id_spath)

# ╔═╡ 7e5c01f5-a1fa-4bca-a7eb-c079a7237a8f
md"to facilitate writing tests, this function maps an atom and bond sequence along a path to the `UInt8` sequence stored to represent the path."

# ╔═╡ a25f2e89-66d6-4e83-8275-0541847af896
function merge(atoms::Array{AtomType}, bonds::Array{BondType})
	@assert length(atoms) == length(bonds) + 1 # for validity
	
	ℓ = length(bonds)
	seq = UInt8.(
		[
			i % 2 == 1 ? atoms[i÷2 + 1] : bonds[i÷2]
			for i in 1:(2*ℓ + 1)
		]
	)

	return max(seq, reverse(seq))
end

# ╔═╡ 7d96254c-eac5-4171-88b4-4f29d7097fc3
md"✅ look at the shortest path between node 8 and node 2."

# ╔═╡ 22711aee-059a-48f4-9f4f-43a89e5809e4
begin
	local spaths = get_shortest_paths(mg, 2, 8)
	local spath = spaths[1]
	
	# it's a unique spath
	@test length(spaths) == 1
	@test spath.w ≈ 1.0
	
	# its length is three.
	@test spath.ℓ == 3

	# its atom-bond label sequence
	@test spath.seq == merge(
		[ATOM_N, ATOM_C, ATOM_C, ATOM_O],
		[BOND_AROMATIC, BOND_AROMATIC, BOND_DOUBLE]
	)

	# check
	@test spath.ρ[1] == 8 && spath.ρ[end] == 2
	
	# 1st, 2nd, last atom have correct label in the sequence
	@test spath.seq[1] == UInt8(ATOM_O)
	@test spath.seq[3] == UInt8(ATOM_C)
	@test spath.seq[end] == UInt8(ATOM_N)
end

# ╔═╡ e376b1ec-9140-421d-8b9b-d9df1d026c4e
md"✅ look at an edge case of ℓ = 0."

# ╔═╡ ef3a311e-7d43-4f48-9c12-1bff16006189
@test get_shortest_paths(mg, 8, 8)[1].seq == UInt8.([ATOM_O])

# ╔═╡ ba616dad-8ac4-4709-8730-cc3127bd5833
md"✅ look at an edge case where there is more than one shorest path."

# ╔═╡ a0c5a12d-6beb-47ec-96f8-ccc091dab7c6
begin
	@test length(get_shortest_paths(mg, 12, 6)) == 2
	@test get_shortest_paths(mg, 12, 6)[1].w ≈ 1/2
end

# ╔═╡ 816c5c89-6fb9-4f75-b8bc-60edda482f33
md"✅ look at shorest path between node 12 and node 13."

# ╔═╡ 8c09e491-c03f-40b7-ab2c-faec377ef630
begin
	local spath = get_shortest_paths(mg, 13, 12)[1]
	@test spath.seq == merge(
		[ATOM_N, ATOM_C, ATOM_O],
		[BOND_AROMATIC, BOND_DOUBLE]
	)
	@test AtomType(spath.seq[1]) == mg.atoms[spath.ρ[1]]
end

# ╔═╡ d6adf42a-d8df-46dd-8357-0b72364a396f
md"
## 🐧 explicit shortest path featurization

for testing the kernel, construct the explicit feature vectors.

💡 the kernel between two molecular graphs should be equal to the dot product of their two explicitly constructed shortest path feature vectors.
"

# ╔═╡ 12b6cb4d-6769-4bd0-9bd1-7dee1a3982e2
begin
	test_mgs = MolGraph.(
		# list of test SMILES here.
		[
			"C1=CN=CN1", 
			"NC(N)CC", 
			"NCCC", 
			"CN1C=NC2=C1C(=O)N(C)C(=O)N2C", 
			"O=C(O)c1ccccc1",
			"CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
		]
	)
	
	for mg in test_mgs
		find_shortest_paths!(mg)
	end
end

# ╔═╡ de3c0011-0559-4380-a5f3-7a58dd664ddc
md"### ordinary shortest path features

💡 each shortest path feature is defined by a length and an un-ordered pair of atom labels.

construct a feature space based on active features among this list of molecules.
"

# ╔═╡ b04958c5-c81c-4150-b192-0a97f96afc5e
fs = FeatureSpace{OrdinarySPF}(test_mgs)

# ╔═╡ 1b95ca2e-14bb-4844-b66c-82f4534076c3
md"e.g. here is the feature vector of one of the molecules."

# ╔═╡ ad8a3f0f-4464-40e7-a3da-ab27e9501fff
build_feature_vector(mg, fs)

# ╔═╡ 1a75285e-281a-47b8-b239-a083fc9bb373
md"✅ test some tricky details to get equality to work."

# ╔═╡ ae0d1430-d5f5-41fa-a92c-179359e81cd7
begin
	f1 = LabelSeqSPF(UInt8.([ATOM_C, BOND_SINGLE, ATOM_C]))
	f2 = LabelSeqSPF(UInt8.([ATOM_C, BOND_SINGLE, ATOM_C]))
	@test f1 == f2
	@test hash(f1) == hash(f2)
end

# ╔═╡ fee56513-d201-4cdc-bb49-d822a60645ce
md"✅ case 1: test the ordinary shortest path features" 

# ╔═╡ 04dc912c-c50f-4845-ba58-c9b845a88cdd
begin
	# construct molecular graph & find it shortest paths
	mg₁ = MolGraph("C1=CN=CN1")
	find_shortest_paths!(mg₁)

	# all shortest paths are unique and thus have weight one
	@test all(sp.w == 1.0 for sp in mg₁.spaths)

	# build ordinary shortest path feature vector
	ϕ₁ = build_feature_vector(mg₁, fs)
	
	# there are 15 active shortest paths
	@test sum(ϕ₁) == 15
	# ℓ = 0: 5
	@test length(filter(sp -> sp.ℓ == 0, mg₁.spaths)) == 5 # atoms
	# ℓ = 1: 5
	@test length(filter(sp -> sp.ℓ == 1, mg₁.spaths)) == 5 # bonds
	# ℓ = 2: 4
	@test length(filter(sp -> sp.ℓ == 2, mg₁.spaths)) == 5
	# ℓ ≥ 3: NONE! cuz SHORTEST PATH is never 3 but 2.
	@test length(filter(sp -> sp.ℓ > 2, mg₁.spaths)) == 0
	
	# the shortest path node 2 <--> node 4
	local sps = get_shortest_paths(mg₁, 2, 4)
	@test length(sps) == 1
	@test sps[1].ρ == reverse([4, 3, 2])
	@test sps[1].ℓ == 2
	
	viz(mg₁, nlabels=true)
end

# ╔═╡ f56d3c19-9b8d-4a6e-bc57-341701d45d91
md"✅ case 2: test the ordinary shortest path features" 

# ╔═╡ 3e9c0f7a-fe4d-4fcc-8e55-b65a4f309393
begin
	mg₂ = MolGraph("NC(N)CC")
	find_shortest_paths!(mg₂)
	ϕ₂ = build_feature_vector(mg₂, fs)
	@test length(filter(sp -> sp.ℓ == 0, mg₂.spaths)) == length(mg₂.atoms) == 5
	@test length(filter(sp -> sp.ℓ == 1, mg₂.spaths)) == length(mg₂.bonds) == 4
 	@test length(filter(sp -> sp.ℓ == 2, mg₂.spaths)) == 4
	@test length(filter(sp -> sp.ℓ == 3, mg₂.spaths)) == 2
	@test sum(ϕ₂) == 5 + 4 + 4 + 2
	
	viz(mg₂, nlabels=true)
end

# ╔═╡ f998cb67-3470-4287-9496-7133bdd2e9e1
md"✅ test molecular graph kernel, manually computed" 

# ╔═╡ 3e630fb9-8fd1-4cce-a0a8-78d60b9abe6a
PlutoUI.LocalResource("manual_spgk_test.png")

# ╔═╡ 357d1312-f98d-4a3e-9672-f540eba61ebd
@test shortest_path_graph_kernel(mg₁, mg₂) == 9 + 4 + 2 + 8 + 2 + 4 + 1

# ╔═╡ 50202a67-095e-47b2-bb96-c42d0f05bad7
# zero edge matches, so just ℓ = 0 paths match
@test shortest_path_graph_kernel(mg₁, mg₂, exact_seq_matching=true) == 9 + 4

# ╔═╡ c83aec87-9f63-408c-839d-86f09d6b45dd
md"✅ the graph kernel is the inner product between the feature vectors!" 

# ╔═╡ 749cc889-8dfc-41c3-b04b-94e671e3d347
@test dot(ϕ₁, ϕ₂) == shortest_path_graph_kernel(mg₁, mg₂)

# ╔═╡ 7b65e9fd-ce90-49e6-9f55-5bfc289dfa23
begin
	mg₃ = MolGraph("NCCC")
	find_shortest_paths!(mg₃)
	
	# ℓ = 0: [CC, NN]. [3, 1] ⋅ [3, 2] = 11
	# ℓ = 1: [CC, CN]. [2, 2] ⋅ [2, 1] = 4 + 2
	# ℓ = 2: [CC, CN, NN]. [1, 2, 1] ⋅ [1, 1, 0] = 3
	# ℓ = 3: [CN]. [1] ⋅ [2] = 2
	@test shortest_path_graph_kernel(mg₂, mg₃) == 11 + 6 + 3 + 2
	
	# now, when label sequence matters, what are mis-matches? none!
	@test shortest_path_graph_kernel(mg₂, mg₃) == shortest_path_graph_kernel(
		mg₂, mg₃, exact_seq_matching=true)
	
	viz(mg₃, nlabels=true)
end

# ╔═╡ 9fe0ea1d-2aa4-42f9-80f0-278e2937cbca
md"✅ a more complex example: caffeine and ibuprofen."

# ╔═╡ 020f859d-0f5a-4b82-bcc5-88dedc1776d2
begin
	caf = MolGraph("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
	ibu = MolGraph("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
	find_shortest_paths!.([caf, ibu])
	
	ϕ_caf = build_feature_vector(caf, fs)
	ϕ_ibu = build_feature_vector(ibu, fs)

	@test dot(ϕ_ibu, ϕ_caf) == shortest_path_graph_kernel(caf, ibu)
end

# ╔═╡ 4c01eaf7-d479-4648-be16-27e7eef7a309
md"✅ for ibuprofen, test some of the shortest paths and shortest path features."

# ╔═╡ 0a0b68fe-c0bb-433f-8cb8-e53465664bd0
begin
	@test get_shortest_paths(ibu, 10, 8)[1].ℓ == 2
	@test get_shortest_paths(ibu, 14, 12)[1].ℓ == 3
	
	# test specific SPFs
	@test ϕ_ibu[
		fs.feature_to_id[OrdinarySPF(ATOM_O, ATOM_O, 2)]
	] == 1
	@test ϕ_ibu[
		fs.feature_to_id[OrdinarySPF(ATOM_C, ATOM_O, 8)]
	] == 2
	@test ϕ_ibu[
		fs.feature_to_id[OrdinarySPF(ATOM_O, ATOM_O, 0)]
	] == 2
	@test ϕ_ibu[
		fs.feature_to_id[OrdinarySPF(ATOM_C, ATOM_O, 1)] 
	]== 2
	
	@test get_shortest_paths(ibu, 12, 9)[1].ρ == [12, 11, 8, 9]
	@test get_shortest_paths(ibu, 5, 3)[1].ρ == reverse([5, 4, 2, 3])

	viz(ibu, nlabels=true, size=(500, 500))
end

# ╔═╡ 6c0bbcce-0d91-4bd8-83e8-6e663eaf10a9
md"✅ test shortest path features for another molecule."

# ╔═╡ 3abe5e9e-777e-44fd-bd3e-8611beacc6f4
begin
	ba = MolGraph("O=C(O)c1ccccc1")
	find_shortest_paths!(ba)
	
	ϕ_ba = build_feature_vector(ba, fs)
	
	# test specific SPFs
	@test ϕ_ba[
		fs.feature_to_id[OrdinarySPF(ATOM_O, ATOM_O, 2)]
	] == 1
	@test ϕ_ba[
		fs.feature_to_id[OrdinarySPF(ATOM_O, ATOM_O, 0)]
	] == 2
	@test ϕ_ba[
		fs.feature_to_id[OrdinarySPF(ATOM_C, ATOM_C, 0)]
	] == 7

	@test get_shortest_paths(ba, 1, 3)[1].ℓ == 2
	@test get_shortest_paths(ba, 2, 8)[1].ℓ == 3
	
	@test length(get_shortest_paths(ba, 2, 7)) == 2

	@test dot(ϕ_ba, ϕ_ibu) == shortest_path_graph_kernel(ba, ibu)
	@test dot(ϕ_ba, ϕ_ba) == shortest_path_graph_kernel(ba, ba)

	@test get_shortest_paths(ba, 1, 3)[1].ρ == [1, 2, 3]
	@test get_shortest_paths(ba, 7, 5)[1].ρ == reverse([7, 6, 5])
	@test get_shortest_paths(ba, 8, 1)[1].ρ == reverse([8, 9, 4, 2, 1])

	# symmetry
	@test shortest_path_graph_kernel(ba, ibu) == shortest_path_graph_kernel(ibu, ba)

	@test shortest_path_graph_kernel(ba, ibu, exact_seq_matching=true) < shortest_path_graph_kernel(ba, ibu, exact_seq_matching=false)
	
	viz(ba, nlabels=true)
end

# ╔═╡ 9f0859f6-7383-4938-919d-4edf1bed5e90
shortest_path_graph_kernel(ba, ibu, exact_seq_matching=true)

# ╔═╡ 1e6fc4f0-a1fe-43e0-8f0b-97f112fdf1c5
md"single ion"

# ╔═╡ 0561aa20-28ef-416a-871e-82b0703efff0
begin
	ion = MolGraph("[Na+]")
	find_shortest_paths!(ion)
	@test length(ion.spaths) == 1

	@test shortest_path_graph_kernel(ion, ion) == 1
end

# ╔═╡ e08d4d33-30fa-4fea-807e-7c8ab244debc
md"### exact matching on the atom-bond label sequence

💡 each feature now is defined by an atom-bond label sequence.
"

# ╔═╡ 8d7e1fc5-0026-4e41-8ea5-8a3dc5877431
fs_seq = FeatureSpace{LabelSeqSPF}(test_mgs)

# ╔═╡ 52da25a4-dd70-4531-aafd-d69ebcef95a7
md"e.g. here is the feature vector of one of the molecules."

# ╔═╡ ad760d7d-c354-43c7-b5c7-c879cf48109b
build_feature_vector(mg, fs_seq)

# ╔═╡ 7d5b9446-9adc-489c-961c-a6b6ea3240ee
begin
	ϕ_seq_ba = build_feature_vector(ba, fs_seq)
	ϕ_seq_ibu = build_feature_vector(ibu, fs_seq)
	
	# just one {O, O} ℓ = 2 feature
	@test ϕ_seq_ba[
		fs_seq.feature_to_id[
			LabelSeqSPF(
				merge(
					[ATOM_O, ATOM_C, ATOM_O],
					[BOND_DOUBLE, BOND_SINGLE]
				)
			)
		]
	] == 1
	@test dot(ϕ_seq_ba, ϕ_seq_ibu) ≈ shortest_path_graph_kernel(
		ba, ibu, exact_seq_matching=true
	)
end

# ╔═╡ 303a7141-a442-4162-b67c-f7d12c3f4ca1
# exact matching is more sparse.
@test shortest_path_graph_kernel(
	ba, ibu, exact_seq_matching=true
) < shortest_path_graph_kernel(
	ba, ibu, exact_seq_matching=false
)

# ╔═╡ cbd9304a-d88e-42a6-8cc9-ad5ef6198f36
begin
	# ethanol vs carbon monoxide
	mg_et = MolGraph("CCO") # [C, SINGLE, C, SINGLE, O]
	mg_co = MolGraph("C#O")
	find_shortest_paths!.([mg_et, mg_co])
	@test shortest_path_graph_kernel(
		mg_et, mg_co, exact_seq_matching=false
	) == 4.0 # 2 x C, 1 x O, 1 x CO

	@test shortest_path_graph_kernel(
		mg_et, mg_co, exact_seq_matching=true
	) == 3.0 # 2 x C, 1 x O, 0 x CO cuz triple vs single blod
end

# ╔═╡ 56f05958-1237-4bcf-a2e4-0092f75bd443
md"### exact vs non-exact matching"

# ╔═╡ d8c7de3d-5452-489c-8228-78ff131ac840
begin
	cnc = MolGraph("CNC")
	find_shortest_paths!(cnc)

	ccc = MolGraph("CCC")
	find_shortest_paths!(ccc)

	# ℓ = 0: 
	#    cnc: C, N, C 
	#    ccc: C, C, C
	# ℓ = 1
	#    cnc: CN, NC
	#    ccc: CC, CC
	# ℓ = 2
	#    cnc: CNC
	#    ccc: CCC
	@test shortest_path_graph_kernel(cnc, ccc, exact_seq_matching=false) == 6 + 0 + 1
	@test shortest_path_graph_kernel(cnc, ccc, exact_seq_matching=true) == 6 + 0 + 0
end

# ╔═╡ 54921516-63b0-472a-b1d3-005a73762a90
md"## 🐧 multi-threaded Gram matrix computation"

# ╔═╡ d2eb0777-0b84-4281-8363-2dbde45cbe0b
Threads.nthreads() # number of threads

# ╔═╡ 70359912-dc86-4c73-ac53-85001f27b65a
K = @btime compute_Gram_matrix(vcat(test_mgs, test_mgs[[3, 1]]), true) # put in duplicates on purpose for later

# ╔═╡ f4c2a738-d82f-445c-890f-cec99fbd3573
@test K[3, 4] ≈ shortest_path_graph_kernel(
	test_mgs[3], test_mgs[4], exact_seq_matching=true
)

# ╔═╡ 3488704a-55d6-40f3-a1b3-86f2feec3fbc
begin
	idps = indistinguishable_pairs(K)
	@test length(idps) == 2
	@test (1, 8) in idps
	@test (3, 7) in idps
	idps
end

# ╔═╡ 27f09dc8-fe7b-4994-b0e5-29de901d5e8c
md"### 🕐 timing"

# ╔═╡ 214f70c0-1a5d-4f6c-a907-8e86025fe19c
@btime shortest_path_graph_kernel(ba, ibu)

# ╔═╡ 2bfb6517-de01-4c4e-bf40-346e18f129c5
@btime shortest_path_graph_kernel(ba, ibu, exact_seq_matching=true)

# ╔═╡ aa200f41-69cf-4384-a54b-079c19b88c21
md"### centering and normalization"

# ╔═╡ 5f31115a-4e7d-4282-a6f4-62d08326772d
begin
	# make fake feature vectors. 10 2D feature vectors
	fXs = rand(2, 10) # feature vectors in cols

	# create centered feature vectors
	fX̂s = fXs .- mean(fXs, dims=2)
	
	# make fake Gram matrix
	fK = fXs' * fXs # not centered
	fK̂ = fX̂s' * fX̂s # centered

	# the centered Gram matrix should match the Gram matrix
	#  computed from scratch, from the centered feature vectors
	@test center_Gram_matrix(fK) ≈ fK̂

	# cosine similarity
	cos_sim = [
		dot(fXs[:, i], fXs[:, j]) / (
			norm(fXs[:, i]) * norm(fXs[:, j])
		) for i = 1:10, j = 1:10
	]

	# cosine similarity matches normalized kernel matrix
	@test normalize_Gram_matrix(fK) ≈ cos_sim

	# remains semi-positive definite
	@test all(x -> x >= -1e-10, eigen(cos_sim).values)
end

# ╔═╡ 9cc543f5-2463-4a1f-906a-8da2b04c92f7
md"## 🎓 for teaching

show two molecular graphs and two paths common to them.
"

# ╔═╡ 0555d22b-33b1-46b6-937b-2dca682476d2
viz(
	caf, 
	id_spath=findfirst(caf.spaths .== get_shortest_paths(caf, 2, 9)),
	nlabels=false
)

# ╔═╡ eb8f8e60-4b9e-40af-849e-d5c7fbc274fe
# adenine 
adn = MolGraph("C1=NC2=NC=NC(=C2N1)N"); find_shortest_paths!(adn)

# ╔═╡ 9df77e29-369d-481d-a2f4-b20ee0854c49
viz(
	adn, 
	id_spath=findfirst(adn.spaths .== get_shortest_paths(adn, 6, 9)),
	nlabels=false
)

# ╔═╡ 3fb5f269-ba3e-496b-abd9-158b59db77da
md"## 🐧 handling non-connected graphs"

# ╔═╡ 329a17cf-e24a-4223-b5a2-d1db4c72bbab
begin
	mg_nc₁ = MolGraph("CC(C(=O)OCC)S(=O)(=O)[O-].[Na+]")
	find_shortest_paths!(mg_nc₁)
	viz(mg_nc₁)
end

# ╔═╡ b0284ce6-d3f5-4bc9-bb05-a9417b51410f
@test length(get_shortest_paths(mg_nc₁, 12)) == 1

# ╔═╡ 2ee527b4-02c7-40ab-8921-d3b514fd7849
begin
	mg_nc₂ = MolGraph("CC(C(=O)OCC)S(=O)(=O)[O]")
	find_shortest_paths!(mg_nc₂)
	viz(mg_nc₂)
end

# ╔═╡ e8d82eab-c449-4a48-93f1-1a505a2f1874
# everything same between 1 and 2 except for Na
@test shortest_path_graph_kernel(mg_nc₁, mg_nc₂) ≈ shortest_path_graph_kernel(mg_nc₂, mg_nc₂)

# ╔═╡ f2205da5-9e55-42a6-98d5-156f03d932e2
# Na-Na score for 2-2 similarity
@test shortest_path_graph_kernel(mg_nc₁, mg_nc₁) ≈ shortest_path_graph_kernel(mg_nc₂, mg_nc₂) + 1

# ╔═╡ Cell order:
# ╠═cecf3058-bb8f-11f0-97f3-bda46249b7c9
# ╟─f9461018-9522-4fe6-a53f-665b3839b6e4
# ╠═50d85bc2-2e1a-424a-84f7-5a22c88759ed
# ╟─f774c2df-1859-443e-9fee-9d490a18c83b
# ╠═e978a83d-f0b3-4dd0-98e6-250124698ec4
# ╠═ffa9d1c9-9534-46fd-aed2-ed549ac46384
# ╠═1f17de1a-5375-4c22-97cd-7dc04645b9d1
# ╟─7e5c01f5-a1fa-4bca-a7eb-c079a7237a8f
# ╠═a25f2e89-66d6-4e83-8275-0541847af896
# ╟─7d96254c-eac5-4171-88b4-4f29d7097fc3
# ╠═22711aee-059a-48f4-9f4f-43a89e5809e4
# ╟─e376b1ec-9140-421d-8b9b-d9df1d026c4e
# ╠═ef3a311e-7d43-4f48-9c12-1bff16006189
# ╟─ba616dad-8ac4-4709-8730-cc3127bd5833
# ╠═a0c5a12d-6beb-47ec-96f8-ccc091dab7c6
# ╟─816c5c89-6fb9-4f75-b8bc-60edda482f33
# ╠═8c09e491-c03f-40b7-ab2c-faec377ef630
# ╟─d6adf42a-d8df-46dd-8357-0b72364a396f
# ╠═12b6cb4d-6769-4bd0-9bd1-7dee1a3982e2
# ╟─de3c0011-0559-4380-a5f3-7a58dd664ddc
# ╠═b04958c5-c81c-4150-b192-0a97f96afc5e
# ╟─1b95ca2e-14bb-4844-b66c-82f4534076c3
# ╠═ad8a3f0f-4464-40e7-a3da-ab27e9501fff
# ╟─1a75285e-281a-47b8-b239-a083fc9bb373
# ╠═ae0d1430-d5f5-41fa-a92c-179359e81cd7
# ╟─fee56513-d201-4cdc-bb49-d822a60645ce
# ╠═04dc912c-c50f-4845-ba58-c9b845a88cdd
# ╟─f56d3c19-9b8d-4a6e-bc57-341701d45d91
# ╠═3e9c0f7a-fe4d-4fcc-8e55-b65a4f309393
# ╟─f998cb67-3470-4287-9496-7133bdd2e9e1
# ╟─3e630fb9-8fd1-4cce-a0a8-78d60b9abe6a
# ╠═357d1312-f98d-4a3e-9672-f540eba61ebd
# ╠═50202a67-095e-47b2-bb96-c42d0f05bad7
# ╟─c83aec87-9f63-408c-839d-86f09d6b45dd
# ╠═749cc889-8dfc-41c3-b04b-94e671e3d347
# ╠═7b65e9fd-ce90-49e6-9f55-5bfc289dfa23
# ╟─9fe0ea1d-2aa4-42f9-80f0-278e2937cbca
# ╠═020f859d-0f5a-4b82-bcc5-88dedc1776d2
# ╟─4c01eaf7-d479-4648-be16-27e7eef7a309
# ╠═0a0b68fe-c0bb-433f-8cb8-e53465664bd0
# ╟─6c0bbcce-0d91-4bd8-83e8-6e663eaf10a9
# ╠═3abe5e9e-777e-44fd-bd3e-8611beacc6f4
# ╠═9f0859f6-7383-4938-919d-4edf1bed5e90
# ╟─1e6fc4f0-a1fe-43e0-8f0b-97f112fdf1c5
# ╠═0561aa20-28ef-416a-871e-82b0703efff0
# ╟─e08d4d33-30fa-4fea-807e-7c8ab244debc
# ╠═8d7e1fc5-0026-4e41-8ea5-8a3dc5877431
# ╟─52da25a4-dd70-4531-aafd-d69ebcef95a7
# ╠═ad760d7d-c354-43c7-b5c7-c879cf48109b
# ╠═7d5b9446-9adc-489c-961c-a6b6ea3240ee
# ╠═303a7141-a442-4162-b67c-f7d12c3f4ca1
# ╠═cbd9304a-d88e-42a6-8cc9-ad5ef6198f36
# ╟─56f05958-1237-4bcf-a2e4-0092f75bd443
# ╠═d8c7de3d-5452-489c-8228-78ff131ac840
# ╟─54921516-63b0-472a-b1d3-005a73762a90
# ╠═d2eb0777-0b84-4281-8363-2dbde45cbe0b
# ╠═70359912-dc86-4c73-ac53-85001f27b65a
# ╠═f4c2a738-d82f-445c-890f-cec99fbd3573
# ╠═3488704a-55d6-40f3-a1b3-86f2feec3fbc
# ╟─27f09dc8-fe7b-4994-b0e5-29de901d5e8c
# ╠═214f70c0-1a5d-4f6c-a907-8e86025fe19c
# ╠═2bfb6517-de01-4c4e-bf40-346e18f129c5
# ╟─aa200f41-69cf-4384-a54b-079c19b88c21
# ╠═5f31115a-4e7d-4282-a6f4-62d08326772d
# ╟─9cc543f5-2463-4a1f-906a-8da2b04c92f7
# ╠═0555d22b-33b1-46b6-937b-2dca682476d2
# ╠═eb8f8e60-4b9e-40af-849e-d5c7fbc274fe
# ╠═9df77e29-369d-481d-a2f4-b20ee0854c49
# ╟─3fb5f269-ba3e-496b-abd9-158b59db77da
# ╠═329a17cf-e24a-4223-b5a2-d1db4c72bbab
# ╠═b0284ce6-d3f5-4bc9-bb05-a9417b51410f
# ╠═2ee527b4-02c7-40ab-8921-d3b514fd7849
# ╠═e8d82eab-c449-4a48-93f1-1a505a2f1874
# ╠═f2205da5-9e55-42a6-98d5-156f03d932e2
