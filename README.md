# SOLPS2IMAS.jl

![Format Check](https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl/actions/workflows/format_check.yml/badge.svg)
![Docs](https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl/actions/workflows/make_docs.yml/badge.svg)
![Tests](https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl/actions/workflows/test.yml/badge.svg)

Utility for loading data from existing SOLPS runs, including custom data in b2time.nc,
from native SOLPS output format into IMAS. For installation and usage instructions, see the [online documentation](https://projecttorreypines.github.io/SOLPS2IMAS.jl/stable). For documentation on under development branch, see [dev online documentation](https://projecttorreypines.github.io/SOLPS2IMAS.jl/dev).

## Installation

SOLPS2IMAS is registered with public repository [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/). For installation:

```
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("SOLPS2IMAS")
```

## Running tests:

Following code ensures it uses the enviroment described in ./Project.toml:

```bash
% julia --project test/runtests.jl
Test Summary:          | Pass  Total  Time
Test index conversions |    1      1  0.0s
Test Summary:              | Pass  Total  Time
Test file parsing in depth |    6      6  0.4s
Test Summary:       | Pass  Total  Time
Test read_b2_output |   19     19  1.9s
solps2imas() time:  49.510687 seconds (154.81 M allocations: 10.512 GiB, 8.27% gc time, 96.47% compilation time)
Test Summary:                        | Pass  Total   Time
Test solps2imas() (overall workflow) |    4      4  50.4s
Test Summary:                                   |  Pass  Total  Time
Test triangular mesh generation from fort files | 23479  23479  3.7s
```
