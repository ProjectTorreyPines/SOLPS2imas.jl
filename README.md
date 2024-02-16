# SOLPS2IMAS.jl

![Format Check](https://github.com/ProjectTorreyPines/SOLPS2IMAS.jl/actions/workflows/format_check.yml/badge.svg)

Utility for loading data from existing SOLPS runs, including custom data in b2time.nc,
from native SOLPS output format into IMAS

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
