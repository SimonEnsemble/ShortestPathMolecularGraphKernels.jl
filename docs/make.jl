using Documenter, ShortestPathMolecularGraphKernels

makedocs(
    sitename="Shortest Path Molecular Graph Kernels",
    repo="github.com/SimonEnsemble/ShortestPathMolecularGraphKernels.jl.git",
    format=Documenter.HTML(; assets=["assets/flux.css"])
)

deploydocs(
    repo="github.com/SimonEnsemble/ShortestPathMolecularGraphKernels.jl.git",
    devbranch="main",
    devurl="stable",   # put docs in stable no matter what to keep main up to date
    push_preview=true  # builds docs for pull requests too
)
