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

function test_read_b2_output()
    filename = "samples/b2fstate"
    contents = read_b2_output(filename)
    nx, ny, ns = contents["nx,ny,ns"]
    nx += 2  # Account for guard cells
    ny += 2  # Account for guard cells
    @assert(size(contents["te"]) == (ny, nx))

    filename = "samples/b2fgmtry"
    contents = read_b2_output(filename)
    nxg, nyg = contents["nx,ny"]
    nxg += 2  # Account for guard cells
    nyg += 2  # Account for guard cells
    sc = size(contents["crx"])
    # println("nx=$(nxg), ny=$(nyg), size(contents['crx'])=$(sc)")
    @assert(size(contents["crx"]) == (4, nyg, nxg))
    @assert(nyg == ny)
    @assert(nxg == nx)

    filename = "samples/b2time.nc"
    contents = read_b2_output(filename)
    nx = contents["nx"]
    ny = contents["ny"]
    nt = contents["time"]
    print(size(contents["ne2d"]))
    print((nt, ny, nx))
    @assert(size(contents["ne2d"]) == (nt, ny, nx))

    return true
end

@testset "omasstuff" begin
    @test try_omas() === nothing
    @test populate_grid_ggd(94, 38) === nothing
    @test test_read_b2_output()
end
