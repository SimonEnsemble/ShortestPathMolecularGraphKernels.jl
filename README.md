# ShortestPathMolecularGraphKernels.jl

[![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://SimonEnsemble.github.io/ShortestPathMolecularGraphKernels.jl/stable)

```julia
caf = MolGraph("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")   # caffeine
adn = MolGraph("C1=NC2=NC=NC(=C2N1)N")           # adenine
find_shortest_paths!.([caf, adn])                 

similarity_score = shortest_path_graph_kernel(caf, adn, exact_seq_matching=false) # 367.0
```

![shortest path molecular graph kernel overview](docs/src/overview.png)
