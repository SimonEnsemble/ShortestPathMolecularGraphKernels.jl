### A Pluto.jl notebook ###
# v0.20.23

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
end

# ╔═╡ f9461018-9522-4fe6-a53f-665b3839b6e4
md"# 👀 demo and test along the way


### building a molecular graph and listing its shortest paths

🍕 construct a molecular graph from the SMILES string, then visualize it.
"

# ╔═╡ 50d85bc2-2e1a-424a-84f7-5a22c88759ed
mg = MolGraph(
	"CN1C=NC2=C1C(=O)N(C)C(=O)N2C", # SMILES
	false,                          # include H atoms?
	ℓ_max=typemax(Int)              # max shortest path length to consider
)

# ╔═╡ f774c2df-1859-443e-9fee-9d490a18c83b
md"id of shortest path $(@bind id_spath PlutoUI.Slider(1:length(mg.spaths)))"

# ╔═╡ 1f17de1a-5375-4c22-97cd-7dc04645b9d1
viz(mg, nlabels=true, id_spath=id_spath, edge_hl_color=:green3)

# ╔═╡ ffa9d1c9-9534-46fd-aed2-ed549ac46384
mg.spaths[id_spath]

# ╔═╡ 7e5c01f5-a1fa-4bca-a7eb-c079a7237a8f
md"🍕 lookit the shortest paths."

# ╔═╡ a25f2e89-66d6-4e83-8275-0541847af896
function intermingle(atom_seq::Array{AtomType}, bond_seq::Array{BondType})
	ℓ = length(bond_seq)
	@assert length(atom_seq) == ℓ + 1
	seq = Vector{UInt8}(undef, 2 * ℓ + 1)
	for i = 1:ℓ+1
		seq[2 * i - 1] = UInt8(atom_seq[i])
		if i ≤ ℓ
			seq[2 * i] = UInt8(bond_seq[i])
		end
	end
	reverse_seq = reverse(seq)
	if isless(seq, reverse_seq)
		return reverse_seq
	else
		return seq
	end
end

# ╔═╡ 22711aee-059a-48f4-9f4f-43a89e5809e4
begin
	# 2 -> 8
	local spaths = get_shortest_paths(mg, 2, 8)
	local spath = spaths[1]
	@test length(spaths) == 1 # unique spath
	@test spath.ℓ == 3
	@test spath.seq == intermingle(
		[ATOM_N, ATOM_C, ATOM_C, ATOM_O],
		[BOND_AROMATIC, BOND_AROMATIC, BOND_DOUBLE]
	)
	@test spath.w ≈ 1.0
	@test spath.ρ[1] == 8 && spath.ρ[end] == 2
	@test spath.seq[1] == UInt8(mg.atoms[spath.ρ[1]])
	@test spath.seq[3] == UInt8(mg.atoms[spath.ρ[2]])
	@test spath.seq[end] == UInt8(mg.atoms[spath.ρ[end]])
	
	# edge case of ℓ = 0
	@test get_shortest_paths(mg, 8, 8)[1].seq == UInt8.([ATOM_O])

	# edge of case of more than one shortest path
	@test length(get_shortest_paths(mg, 12, 6)) == 2
	@test get_shortest_paths(mg, 12, 6)[1].w ≈ 1/2

	# 13 -> 12
	local spath = get_shortest_paths(mg, 13, 12)[1]
	@test spath.seq == intermingle(
		[ATOM_N, ATOM_C, ATOM_O],
		[BOND_AROMATIC, BOND_DOUBLE]
	)
	#@test 
	@test AtomType(spath.seq[1]) == mg.atoms[spath.ρ[1]]
end

# ╔═╡ d6adf42a-d8df-46dd-8357-0b72364a396f
md"
### shortest path featurization

🍕 construct the explicit feature vectors"

# ╔═╡ 12b6cb4d-6769-4bd0-9bd1-7dee1a3982e2
test_mgs = MolGraph.(
	# list of test SMILES here.
	[
		"C1=CN=CN1", 
		"NC(N)CC", 
		"NCCC", 
		"CN1C=NC2=C1C(=O)N(C)C(=O)N2C", 
		"O=C(O)c1ccccc1",
		"CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
		"CN1C=NC2=C1C(=O)N(C)C(=O)N2C"
	], 
	false
)

# ╔═╡ b04958c5-c81c-4150-b192-0a97f96afc5e
fs = FeatureSpace{OrdinarySPF}(test_mgs)

# ╔═╡ ad8a3f0f-4464-40e7-a3da-ab27e9501fff
build_feature_vector(mg, fs)

# ╔═╡ 8d7e1fc5-0026-4e41-8ea5-8a3dc5877431
fs_seq = FeatureSpace{LabelSeqSPF}(test_mgs)

# ╔═╡ ad760d7d-c354-43c7-b5c7-c879cf48109b
build_feature_vector(mg, fs_seq)

# ╔═╡ 1a75285e-281a-47b8-b239-a083fc9bb373
md"tricky details to get sets to work."

# ╔═╡ ae0d1430-d5f5-41fa-a92c-179359e81cd7
begin
	f1 = LabelSeqSPF(UInt8.([ATOM_C, BOND_SINGLE, ATOM_C]))
	f2 = LabelSeqSPF(UInt8.([ATOM_C, BOND_SINGLE, ATOM_C]))
	@test f1 == f2
	@test hash(f1) == hash(f2)
end

# ╔═╡ f32fc3fb-c7cb-4154-af41-fc690d9ad28c
md"# tests

🍕 an example I sketched out in my notebook.
"

# ╔═╡ 04dc912c-c50f-4845-ba58-c9b845a88cdd
begin
	mg₁ = MolGraph("C1=CN=CN1", false)
	ϕ₁ = build_feature_vector(mg₁, fs)
	@test sum(ϕ₁) == 15
	local sps = get_shortest_paths(mg₁, 2, 4)
	@test length(sps) == 1
	@test sps[1].ρ == reverse([4, 3, 2])
	@test sps[1].ℓ == 2
	@test all(sp.w == 1.0 for sp in mg₁.spaths)
	viz(mg₁, nlabels=true)
end

# ╔═╡ 3e9c0f7a-fe4d-4fcc-8e55-b65a4f309393
begin
	mg₂ = MolGraph("NC(N)CC", false)
	ϕ₂ = build_feature_vector(mg₂, fs)
	@test sum(ϕ₂) == 15
	viz(mg₂, nlabels=true)
end

# ╔═╡ 357d1312-f98d-4a3e-9672-f540eba61ebd
@test shortest_path_graph_kernel(mg₁, mg₂) == 9 + 4 + 2 + 8 + 2 + 4 + 1

# ╔═╡ 50202a67-095e-47b2-bb96-c42d0f05bad7
@test shortest_path_graph_kernel(mg₁, mg₂, exact_seq_matching=true) == 9 + 4 # no edge matches, so...

# ╔═╡ 749cc889-8dfc-41c3-b04b-94e671e3d347
@test dot(ϕ₁, ϕ₂) == shortest_path_graph_kernel(mg₁, mg₂)

# ╔═╡ 7b65e9fd-ce90-49e6-9f55-5bfc289dfa23
begin
	mg₃ = MolGraph("NCCC", false)
	# reasoning:
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
md"🍕 caffeine and ibuprofen---compare with featurization."

# ╔═╡ 020f859d-0f5a-4b82-bcc5-88dedc1776d2
begin
	caf = MolGraph("CN1C=NC2=C1C(=O)N(C)C(=O)N2C", false)
	ibu = MolGraph("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", false)
	
	ϕ_caf = build_feature_vector(caf, fs)
	ϕ_ibu = build_feature_vector(ibu, fs)

	@test dot(ϕ_ibu, ϕ_caf) == shortest_path_graph_kernel(caf, ibu)
end

# ╔═╡ 64d79622-dff2-4e77-85d2-4aa782144924
active_features([ibu], OrdinarySPF)

# ╔═╡ 4c01eaf7-d479-4648-be16-27e7eef7a309
md"🍕 for ibuprofen, test some of the shortest paths."

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
end

# ╔═╡ 3c0113a0-6453-4d5b-9d68-614e975d0ece
viz(ibu, nlabels=true)

# ╔═╡ 5d180fdf-e1c0-44f0-8edd-573c9afb7475
viz(caf, nlabels=true)

# ╔═╡ cf8998e1-f39a-45e2-b7e5-1f831cb3a586
md"🍕 couple more..."

# ╔═╡ 3abe5e9e-777e-44fd-bd3e-8611beacc6f4
begin
	ba = MolGraph("O=C(O)c1ccccc1", false)
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
	@test length(get_shortest_paths(ba, 2, 7)) == 2
	viz(ba, nlabels=true)
end

# ╔═╡ a85711a1-8dbd-4bf8-9046-1cdf74a88cd5
begin
	# just one {O, O} ℓ = 2 feature
	@test ϕ_ba[
		fs.feature_to_id[OrdinarySPF(ATOM_O, ATOM_O, 2)]
	] == 1
	
	@test get_shortest_paths(ba, 1, 3)[1].ℓ == 2
	@test get_shortest_paths(ba, 2, 8)[1].ℓ == 3

	@test dot(ϕ_ba, ϕ_ibu) == shortest_path_graph_kernel(ba, ibu)
	# symmetry
	@test shortest_path_graph_kernel(ba, ibu) == shortest_path_graph_kernel(ibu, ba)
	@test dot(ϕ_ba, ϕ_ba) == shortest_path_graph_kernel(ba, ba)

	@test get_shortest_paths(ba, 1, 3)[1].ρ == [1, 2, 3]
	@test get_shortest_paths(ba, 7, 5)[1].ρ == reverse([7, 6, 5])
	@test get_shortest_paths(ba, 8, 1)[1].ρ == reverse([8, 9, 4, 2, 1])
end

# ╔═╡ 88946927-f014-4bec-a6a5-17ed9db02b19
viz(ibu, nlabels=true)

# ╔═╡ 57215f37-7b88-4edf-9411-91cd696238fd
shortest_path_graph_kernel(ba, ba)

# ╔═╡ 214f70c0-1a5d-4f6c-a907-8e86025fe19c
@btime shortest_path_graph_kernel(ba, ibu)

# ╔═╡ 2bfb6517-de01-4c4e-bf40-346e18f129c5
shortest_path_graph_kernel(ba, ba)

# ╔═╡ efdeeed8-9583-4433-b9e0-89aa1ef888d8
@test shortest_path_graph_kernel(ba, ibu, exact_seq_matching=true) < shortest_path_graph_kernel(ba, ibu, exact_seq_matching=false)

# ╔═╡ e08d4d33-30fa-4fea-807e-7c8ab244debc
md"## exact pattern matching

first, get a list of the unique features.
"

# ╔═╡ 7d5b9446-9adc-489c-961c-a6b6ea3240ee
begin
	ϕ_seq_ba = build_feature_vector(ba, fs_seq)
	ϕ_seq_ibu = build_feature_vector(ibu, fs_seq)
	# just one {O, O} ℓ = 2 feature
	@test ϕ_seq_ba[
		fs_seq.feature_to_id[
			LabelSeqSPF(
				intermingle(
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
	mg_et = MolGraph("CCO", false) # [C, SINGLE, C, SINGLE, O]
	mg_co = MolGraph("C#O", false)
	@test shortest_path_graph_kernel(
		mg_et, mg_co, exact_seq_matching=false
	) == 4.0 # 2 x C, 1 x O, 1 x CO

	@test shortest_path_graph_kernel(
		mg_et, mg_co, exact_seq_matching=true
	) == 3.0 # 2 x C, 1 x O, 0 x CO cuz triple vs single blod
end

# ╔═╡ Cell order:
# ╠═cecf3058-bb8f-11f0-97f3-bda46249b7c9
# ╟─f9461018-9522-4fe6-a53f-665b3839b6e4
# ╠═50d85bc2-2e1a-424a-84f7-5a22c88759ed
# ╟─f774c2df-1859-443e-9fee-9d490a18c83b
# ╠═1f17de1a-5375-4c22-97cd-7dc04645b9d1
# ╠═ffa9d1c9-9534-46fd-aed2-ed549ac46384
# ╟─7e5c01f5-a1fa-4bca-a7eb-c079a7237a8f
# ╠═a25f2e89-66d6-4e83-8275-0541847af896
# ╠═22711aee-059a-48f4-9f4f-43a89e5809e4
# ╟─d6adf42a-d8df-46dd-8357-0b72364a396f
# ╠═12b6cb4d-6769-4bd0-9bd1-7dee1a3982e2
# ╠═b04958c5-c81c-4150-b192-0a97f96afc5e
# ╠═ad8a3f0f-4464-40e7-a3da-ab27e9501fff
# ╠═64d79622-dff2-4e77-85d2-4aa782144924
# ╠═8d7e1fc5-0026-4e41-8ea5-8a3dc5877431
# ╠═ad760d7d-c354-43c7-b5c7-c879cf48109b
# ╟─1a75285e-281a-47b8-b239-a083fc9bb373
# ╠═ae0d1430-d5f5-41fa-a92c-179359e81cd7
# ╟─f32fc3fb-c7cb-4154-af41-fc690d9ad28c
# ╠═04dc912c-c50f-4845-ba58-c9b845a88cdd
# ╠═3e9c0f7a-fe4d-4fcc-8e55-b65a4f309393
# ╠═357d1312-f98d-4a3e-9672-f540eba61ebd
# ╠═50202a67-095e-47b2-bb96-c42d0f05bad7
# ╠═749cc889-8dfc-41c3-b04b-94e671e3d347
# ╠═7b65e9fd-ce90-49e6-9f55-5bfc289dfa23
# ╟─9fe0ea1d-2aa4-42f9-80f0-278e2937cbca
# ╠═020f859d-0f5a-4b82-bcc5-88dedc1776d2
# ╟─4c01eaf7-d479-4648-be16-27e7eef7a309
# ╠═0a0b68fe-c0bb-433f-8cb8-e53465664bd0
# ╠═3c0113a0-6453-4d5b-9d68-614e975d0ece
# ╠═5d180fdf-e1c0-44f0-8edd-573c9afb7475
# ╟─cf8998e1-f39a-45e2-b7e5-1f831cb3a586
# ╠═3abe5e9e-777e-44fd-bd3e-8611beacc6f4
# ╠═a85711a1-8dbd-4bf8-9046-1cdf74a88cd5
# ╠═88946927-f014-4bec-a6a5-17ed9db02b19
# ╠═57215f37-7b88-4edf-9411-91cd696238fd
# ╠═214f70c0-1a5d-4f6c-a907-8e86025fe19c
# ╠═2bfb6517-de01-4c4e-bf40-346e18f129c5
# ╠═efdeeed8-9583-4433-b9e0-89aa1ef888d8
# ╟─e08d4d33-30fa-4fea-807e-7c8ab244debc
# ╠═7d5b9446-9adc-489c-961c-a6b6ea3240ee
# ╠═303a7141-a442-4162-b67c-f7d12c3f4ca1
# ╠═cbd9304a-d88e-42a6-8cc9-ad5ef6198f36
