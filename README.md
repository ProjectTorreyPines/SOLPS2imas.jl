# SOLPS2IMAS.jl

Utility for loading data from existing SOLPS runs, including custom data in b2time.nc, from native SOLPS output format into IMAS

## To-do:

* Test using new [h5i2imas functionality](https://github.com/ProjectTorreyPines/OMAS.jl/commit/904562b4040c857260d832747d1ba46a2bf614f6) added to OMAS.
* Test for more groups and quantities in populate_grid_ggd function. Identify the fields and list them here.
* Expand testing to some self-consistency test of read data other than just running the function without errors.

## Running tests:

Following code ensures it uses the enviroment described in ./Project.toml:

```
(base) gupta@F-CJXNMY7L7 SOLPS2IMAS.jl % julia --project=./ ./test/runtests.jl 
it's the omas function
another fun function!!!!!11!!!!!
Test Summary: | Pass  Total  Time
omasstuff     |    3      3  2.5s
```
