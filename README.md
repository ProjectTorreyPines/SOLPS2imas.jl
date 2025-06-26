# SOLPS2imas.jl

![Format Check](https://github.com/ProjectTorreyPines/SOLPS2imas.jl/actions/workflows/format_check.yml/badge.svg)
![Docs](https://github.com/ProjectTorreyPines/SOLPS2imas.jl/actions/workflows/make_docs.yml/badge.svg)
![Tests](https://github.com/ProjectTorreyPines/SOLPS2imas.jl/actions/workflows/test.yml/badge.svg)
[![codecov](https://codecov.io/gh/ProjectTorreyPines/SOLPS2imas.jl/graph/badge.svg?token=ZJBRLAXIS1)](https://codecov.io/gh/ProjectTorreyPines/SOLPS2imas.jl)

Utility for loading data from existing SOLPS runs, including custom data in b2time.nc,
from native SOLPS output format into IMAS. For installation and usage instructions, see the [online documentation](https://projecttorreypines.github.io/SOLPS2imas.jl/stable). For documentation on under development branch, see [dev online documentation](https://projecttorreypines.github.io/SOLPS2imas.jl/dev).

## Installation

```
using Pkg
Pkg.add("SOLPS2imas")
```

## For developers

If you are contributing to this project, you would need to install [dvc](https://dvc.org/) to fetch sample files for testing. Once installed, please configure your ssh so that you can ssh into omega tunneling through cybele without requiring to enter password. This is optional but will make it much easier for you.

Once you have completed above steps, inside the git repo, simply do:
```bash
dvc pull
```

This would download the sample files in the `samples` directory. Then to run tests, you would first need to instantiate the project:
```bash
julia --project
```
Then press `]`:
```julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.11.5 (2025-04-14)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> ]
```
Then type:
```julia
(SOLPS2imas) pkg> instantiate
```
Once the package has been instantiated, you can run the tests using:
```julia
(SOLPS2imas) pkg> test
```
This would run all the tests though. To run specific tests, you can do following from the command line to see help options (this works after you ahve instantiated the project like mentioned above):
```bash
% julia --project test/test.jl help
Usage (from inside SOLPS2imas.jl): 
julia --project test/test.jl [ind] [b2] [solps2imas] [parser] [fort] [boundary] [h] [help]

Run tests. Default is all tests.

Optional arguments:
    ind        Test index conversions                               
    b2         Test pread_b2_output()                               
    solps2imas Test solps2imas                                      
    parser     Stress test file parsing (other than b2 output files)
    fort       Test triangular mesh generation from fort files      
    boundary   Test parsing of boundary parameters                  
    h          Show this help message and exit                      
    help       Show this help message and exit                      

```
