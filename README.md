# SOLPS2IMAS.jl

Utility for loading data from existing SOLPS runs, including custom data in b2time.nc, from native SOLPS output format into IMAS

## OMAS structure notes

1. edge_profiles.ggd store data of quantities; edge_profiles.grid_ggd store data about the grid geometry
2. For both ggd and grid_ggd, different slices are different slices in time. That is, what we were calling grid_number is actually time iteration number/index (say it). ggd[it].time and grid_ggd[it].time will store the time array value.
3. grid_ggd[it]
   1. grid_ggd[it].index needs to be selected based on grid type (1 for linear, 2 for cylinder, 4 for single null etc.)
   2. grid_ggd[it].grid_subset[xx] are various subsets like nodes, edges, cells etc. each with required index value at  grid_ggd[it].grid_subset[xx].index (using xx for indices that do not matter)
   3. grid_ggd[it].grid_subset[xx] properties are automatically generated for nodes, edges or cells but are required to be given for other kinds of subsets (like boundaries)
4. ggd[it]
   1. ggd[it].electrons.<quantity>[xx].grid_subset_index needs to be mentioned to indicate whether the quantity is at a node, edge, or cell.
   2. ggd[it].electrons.<quantity>[xx].values store the actual value.
   4. ggd[it].electrons.<quantity>[xx].coefficients store interpolation function coefficients
   5. For multi-specie data, ggd[it].ions[xx].<quantity>[xx].values is used. For each species, ggd[it].ions[xx].element[:] is used to describe Z and charge state.

## To-do:

- [ ] Read atleast all electron properties.
- [ ] Final function should accept one state/time file and one geometry file to return a OMAS.dd object.
- [ ] Add support for boundary properties.
- [ ] Add support for multi-species data.
- [ ] Test using new [h5i2imas functionality](https://github.com/ProjectTorreyPines/OMAS.jl/commit/904562b4040c857260d832747d1ba46a2bf614f6) added to OMAS.
- [ ] Expand testing to some self-consistency test of read data.

## Running tests:

Following code ensures it uses the enviroment described in ./Project.toml:

```
(base) gupta@F-CJXNMY7L7 SOLPS2IMAS.jl % julia --project=./ ./test/runtests.jl 
it's the omas function
another fun function!!!!!11!!!!!
Test Summary: | Pass  Total  Time
omasstuff     |    3      3  2.5s
```
