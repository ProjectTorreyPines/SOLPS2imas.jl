using SOLPS2IMAS
using Test
using Plots

#To run this test
#$ export LD_LIBRARY_PATH=""
#$ julia
#] add Revise
#] activate .
#] test
# OR
# julia> include("test/runtests.jl")

function test_generate_test_data()
    cr, cz = generate_test_data()
    plot(cr, cz)
    return true
end

@testset "omasstuff" begin
    @test try_omas() === nothing
    @test populate_grid_ggd(94, 38) === nothing
    @test test_generate_test_data()
end
