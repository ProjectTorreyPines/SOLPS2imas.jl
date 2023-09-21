# SOLPS2IMAS.jl

![Format Check](https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl/actions/workflows/format_check.yml/badge.svg)

Utility for loading data from existing SOLPS runs, including custom data in b2time.nc,
from native SOLPS output format into IMAS

## Running tests:

Following code ensures it uses the enviroment described in ./Project.toml:

```
(base) gupta@F-CJXNMY7L7 SOLPS2IMAS.jl % julia --proje
ct=./ test/runtests.jl 
solps2imas() time:   3.309287 seconds (9.92 M allocations: 666.357 MiB, 5.68% gc time, 71.71% compilation time)
Test Summary: | Pass  Total  Time
omasstuff     |    4      4  6.9s
hello there
Test Summary: | Pass  Total  Time
utilities     |    3      3  0.0s
```
