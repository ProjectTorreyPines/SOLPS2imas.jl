# SOLPS2IMAS.jl

![Format Check](https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl/actions/workflows/format_check.yml/badge.svg)
![Docs](https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl/actions/workflows/make_docs.yml/badge.svg)
![Tests](https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl/actions/workflows/test.yml/badge.svg)

Utility for loading data from existing SOLPS runs, including custom data in b2time.nc,
from native SOLPS output format into IMAS. For installation and usage instructions, see the [online documentation](https://projecttorreypines.github.io/SOLPS2IMAS.jl/stable). For documentation on under development branch, see [dev online documentation](https://projecttorreypines.github.io/SOLPS2IMAS.jl/dev).

## Building julia environment (installation)

After cloning this repo, check the make menu:
```
SOLPS2IMAS.jl % make help
Help Menu

make env_with_cloned_repo (or make r): Creates a Julia environment with the cloned repositories
make env_with_git_url (or make u): Creates a Julia environment with the git urls without creating local clones
make clean: Deletes Project.toml and Manifest.toml for a fresh start
```

### make r
This option creates local copies of required private repositories at the same level as current repository and uses them in develop mode to create a Manifest.toml

### make u
This option uses url of required private repositories to create a static Manifest.toml attached to current master branches of these repositories.

### make clean
Deletes Manifest.toml so that environment can be recreated, to update or change the last used method.

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
