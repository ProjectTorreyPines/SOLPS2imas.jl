
# SOLPS2IMAS.jl 

```@contents
Pages = ["index.md"]
Depth = 5
```

## Installation

### Using make:
After cloning this repo, check the make menu:
```
SOLPS2IMAS.jl % make help
Help Menu

make env_with_cloned_repo (or make r): Creates a Julia environment with the cloned repositories
make env_with_git_url (or make u): Creates a Julia environment with the git urls without creating local clones
make clean: Deletes Project.toml and Manifest.toml for a fresh start
```

#### make r
This option creates local copies of required private repositories at the same level as current repository and uses them in develop mode to create a Manifest.toml

#### make u
This option uses url of required private repositories to create a static Manifest.toml attached to current master branches of these repositories.

#### make clean
Deletes Manifest.toml so that environment can be recreated, to update or change the last used method.

### Using Julia REPL and installing using Github url

Or, in julia REPL:
```julia
julia> using Pkg;
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/IMASDD.jl.git");
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl.git");
julia> Pkg.instantiate()
```

## solps2imas

The main function of this module is solps2imas which reads solps input and output files and creates an IMAS IDS instance.

```@docs
solps2imas
```

## Parsing tools

This module uses some parsing functions which can be used standalone as well
```@docs
read_b2_output
read_b2mn_output
read_b2time_output
```

## SOLPS Mesh Tools

```@docs
SOLPS2IMAS.xytoc
SOLPS2IMAS.ctoxy
SOLPS2IMAS.data_xytoc
SOLPS2IMAS.search_point
SOLPS2IMAS.search_edges
SOLPS2IMAS.distance_between_nodes
SOLPS2IMAS.neighbour_inds
SOLPS2IMAS.get_neighbour_inds
SOLPS2IMAS.attach_neightbours!
```

## Subset Identification Tools

These functions provide a way to identify if a linear cell index in SOLPS notation belongs to a grid subset.

```@docs
SOLPS2IMAS.in_core
SOLPS2IMAS.in_sol
SOLPS2IMAS.in_idr
SOLPS2IMAS.in_odr
SOLPS2IMAS.is_x_aligned
SOLPS2IMAS.is_y_aligned
SOLPS2IMAS.is_core_cut
SOLPS2IMAS.is_pfr_cut
SOLPS2IMAS.is_outer_throat
SOLPS2IMAS.is_inner_throat
SOLPS2IMAS.is_outer_midplane
SOLPS2IMAS.is_inner_midplane
SOLPS2IMAS.is_outer_target
SOLPS2IMAS.is_inner_target
SOLPS2IMAS.is_core_boundary
SOLPS2IMAS.is_separatrix
SOLPS2IMAS.get_xpoint_nodes
```