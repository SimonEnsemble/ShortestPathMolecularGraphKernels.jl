"""
    viz(
        mg::MolGraph; id_spath::Union{Int, Nothing}=nothing,
        nlabels::Bool=true, size::Tuple{Int, Int}=(375, 375)
    )

visualize a molecular graph and, if the index of a shortest path is provided,
highlight that shortest path.

uses `GraphMakie.jl` under the hood and displays a plot in a `Pluto.jl` notebook.

# example
```julia
mg = MolGraph("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
viz(mg, id_spath=39) # displays plot in Pluto
```
"""
function viz(
	mg::MolGraph;
	nlabels::Bool=true,
    node_size::Int=25,
    id_spath::Union{Int, Nothing}=nothing,
    edge_hl_color::Symbol=:green3,
    size::Tuple{Int, Int}=(375, 375)
)
    # error out if shortest paths not computed yet
    if length(mg.spaths) == 0 && ! isnothing(id_spath)
        error("$(mg.smiles) lacks shortest paths assigned.")  
    end

	# compute graph layout (re-compute here to avoid storing extra data)
	mol = smilestomol(mg.smiles)
    if mg.incl_hydrogen
        add_hydrogens!(mol)
    end
	coordgen!(mol)
	coords = coords2d(mol).coords

    # for highlighting a path
    default_ecolor = GraphMakie.default_theme(nothing, GraphMakie.LineSegments).color
    edge_color = [default_ecolor for i in 1:ne(mg.g)]
    if ! isnothing(id_spath)
        if mg.spaths[id_spath].ℓ > 0
            for i = 1:mg.spaths[id_spath].ℓ # hops along the path
                u = mg.spaths[id_spath].ρ[i]
                v = mg.spaths[id_spath].ρ[i+1]
                id_ed = get_edge_id(mg.g, u, v)
                edge_color[id_ed] = edge_hl_color
            end
        end
    end
	
	fig = Figure(size=size)
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
		edge_width=3,
		nlabels_distance=12.0,
        edge_color=edge_color
	)
    # increase padding
    autolimits!(ax)
    lims = ax.finallimits[]
    ox, oy = lims.origin
    wx, wy = lims.widths
    pad = 0.3
    limits!(ax, ox-pad, ox+wx+pad, oy-pad, oy+wy+pad)
    #limits!(ax, minimum(lims) - 0.3, maximum(lims) + 0.3)
	
	fig
end
