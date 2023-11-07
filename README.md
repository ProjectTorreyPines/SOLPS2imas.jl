# SOLPS2IMAS.jl

![Format Check](https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl/actions/workflows/format_check.yml/badge.svg)

Utility for loading data from existing SOLPS runs, including custom data in b2time.nc,
from native SOLPS output format into IMAS

## Running tests:

Following code ensures it uses the enviroment described in ./Project.toml:

```
% julia --project test/runtests.jl
Test Summary:          | Pass  Total  Time
Test index conversions |    1      1  0.1s
Test Summary:       | Pass  Total  Time
Test read_b2_output |   19     19  1.8s
solps2imas() time:   8.995446 seconds (16.05 M allocations: 1.115 GiB, 4.49% gc time, 83.09% compilation time: <1% of which was recompilation)
Test Summary:     | Pass  Total  Time
Test solps2imas() |    4      4  9.9s
```
