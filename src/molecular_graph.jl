struct MolGraph
    smiles::String
    incl_hydrogen::Bool
    g::SimpleGraph
    atoms::Vector{AtomType}
    bonds::Vector{BondType}
    spaths::Vector{ShortestPath}
end

function MolGraph(smiles::String, incl_hydrogen::Bool; verbose::Bool=false, ℓ_max::Int=typemax(Int))
    if verbose
        println("reading in ", smiles)
    end

    # use MolecularGraph.jl to read SMILES
    mol = smilestomol(smiles)
    if incl_hydrogen
        add_hydrogens!(mol)
    end

    # construct simple graph
    g = deepcopy(mol.graph)
    if ! is_connected(g)
        error(
            "graph for $smiles is not connected!
            filter ur molecules via:
            [is_connected(smilestomol(smiles).graph) for smiles in my_smiles]
            "
        )
    end
    
    # node labels: atom type
    atoms = [AtomType(elements[mol.vprops[v].symbol].number) for v = 1:nv(g)]
    
    # edge labels: bond order
    aromatic_bond = is_edge_aromatic(mol)
    bonds = [
        aromatic_bond[i] ? BondType(0) : BondType(mol.eprops[ed].order)
            for (i, ed) in enumerate(edges(g))
    ]

    # shortest path lengths between all pairs of graphs
    spaths = find_all_shortest_paths(g, atoms, bonds, ℓ_max=ℓ_max)
    
    return MolGraph(smiles, incl_hydrogen, g, atoms, bonds, spaths)
end

function viz(
	mg::MolGraph;
	nlabels::Bool=false,
    node_size::Int=25,
    id_spath::Union{Int, Nothing}=nothing,
    edge_hl_color::Symbol=:red
)
	# compute layout (re-compute here to avoid storing extra data)
	mol = smilestomol(mg.smiles)
    if mg.incl_hydrogen
        add_hydrogens!(mol)
    end
	coordgen!(mol)
	coords = coords2d(mol).coords

    # edge colors
    default_ecolor = GraphMakie.default_theme(nothing, GraphMakie.LineSegments).color
    edge_color = [default_ecolor for i in 1:ne(mg.g)]
    if ! isnothing(id_spath)
        if mg.spaths[id_spath].ℓ > 0
            for i = 1:mg.spaths[id_spath].ℓ
                id_ed = get_edge_id(mg.g, mg.spaths[id_spath].ρ[i], mg.spaths[id_spath].ρ[i+1])
                edge_color[id_ed] = edge_hl_color
            end
        end
    end
	
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=DataAspect())
	hidedecorations!(ax)
	hidespines!(ax)
	graphplot!(
		mg.g, 
		ilabels=["$(elements[Int(a)].symbol)" for a in mg.atoms],
		nlabels=nlabels ? ["$i" for i = 1:nv(mg.g)] : nothing,
		elabels=["$(Int(o))" for o in mg.bonds],
		node_size=node_size,
		node_color=[elements[Int(a)].cpk_hex for a in mg.atoms],
		layout=coords,
		node_strokewidth=2,
		edge_width=2,
		nlabels_distance=12.0,
        edge_color=edge_color
	)
	
	fig
end

get_shortest_paths(mg::MolGraph, u::Int, v::Int) = get_shortest_paths(mg.spaths, u, v)
