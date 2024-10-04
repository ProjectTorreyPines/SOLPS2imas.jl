
# SOLPS2imas.jl 

```@contents
Pages = ["index.md"]
Depth = 5
```

## Installation

For installation:

```
using Pkg
Pkg.add("SOLPS2imas")
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
```

## SOLPS Mesh Tools

```@docs
SOLPS2imas.xytoc
SOLPS2imas.ctoxy
SOLPS2imas.data_xytoc
SOLPS2imas.search_point
SOLPS2imas.search_edges
SOLPS2imas.distance_between_nodes
SOLPS2imas.neighbour_inds
SOLPS2imas.get_neighbour_inds
SOLPS2imas.attach_neightbours!
```

## Subset Identification Tools

These functions provide a way to identify if a linear cell index in SOLPS notation belongs to a grid subset.

```@docs
SOLPS2imas.in_core
SOLPS2imas.in_sol
SOLPS2imas.in_idr
SOLPS2imas.in_odr
SOLPS2imas.is_x_aligned
SOLPS2imas.is_y_aligned
SOLPS2imas.is_core_cut
SOLPS2imas.is_pfr_cut
SOLPS2imas.is_outer_throat
SOLPS2imas.is_inner_throat
SOLPS2imas.is_outer_midplane
SOLPS2imas.is_inner_midplane
SOLPS2imas.is_outer_target
SOLPS2imas.is_inner_target
SOLPS2imas.is_core_boundary
SOLPS2imas.is_separatrix
SOLPS2imas.get_xpoint_nodes
```