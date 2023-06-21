using SOLPS2IMAS
using Test

#To run this test
#$ export LD_LIBRARY_PATH=""
#$ julia
#] add Revise
#] activate .
#] test
# OR
# julia> include("test/runtests.jl")

@testset "omasstuff" begin
    @test try_omas() === nothing
    @test populate_grid_ggd() === nothing
end
