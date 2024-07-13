
# SOLPS2IMAS.jl 

```@contents
Pages = ["index.md"]
Depth = 5
```

## Installation

SOLPS2IMAS is registered with public repository [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/). For installation:

```
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("SOLPS2IMAS")
```
```

## solps2imas

The main function of this module is solps2imas which reads SOLPS input and output files and creates an IMAS IDS instance.

```@docs
solps2imas
```

The main solps2imas feature can be called to load everything. However, some special cases may demand only subsets of the data. In these cases, it is possible to call subcomponents of the SOLPS to IMAS data transfer feature by using these functions:

```@docs
load_summary_data!
```

## Parsing tools

This module uses some parsing functions which can be used standalone as well
```@docs
read_b2_output
read_b2mn_output
read_b2time_output
read_b2_boundary_parameters
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