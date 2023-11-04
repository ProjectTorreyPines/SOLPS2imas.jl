using Pkg
Pkg.activate(".")
Pkg.rm("OMAS")
Pkg.add(; url="git@github.com:ProjectTorreyPines/OMAS.jl.git")
Pkg.instantiate()