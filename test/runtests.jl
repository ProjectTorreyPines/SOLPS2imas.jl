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
    @assert(size(contents["na"]) == (ns, ny, nx))
    @assert(size(contents["fna"]) == (ns, 2, ny, nx))
    @assert(size(contents["fhe"]) == (ns, ny, nx))

    filename = "samples/b2fgmtry"
    contents = read_b2_output(filename)
    nxg, nyg = contents["nx,ny"]
    nxg += 2  # Account for guard cells
    nyg += 2  # Account for guard cells
    sc = size(contents["crx"])
    @assert(size(contents["crx"]) == (4, nyg, nxg))
    @assert(nyg == ny)
    @assert(nxg == nx)

    filename = "samples/b2time.nc"
    contents = read_b2_output(filename)
    nx = contents["nx"]
    ny = contents["ny"]
    nybl = contents["nybl"]
    nybr = contents["nybr"]
    nt = contents["time"]
    ns = contents["ns"]
    @assert(size(contents["ne2d"]) == (nt, ny, nx))
    @assert(size(contents["ft3dl"]) == (nt, nybl))
    @assert(size(contents["fl3dl"]) == (nt, nybl))
    @assert(size(contents["fc3dl"]) == (nt, nybl))
    @assert(size(contents["ft3dr"]) == (nt, nybr))
    @assert(size(contents["fl3dr"]) == (nt, nybr))
    @assert(size(contents["fc3dr"]) == (nt, nybr))
    @assert(size(contents["fna3da"]) == (nt, ns, nybr))
    @assert(size(contents["tmhacore"]) == (nt, ))
    @assert(size(contents["tmhasol"]) == (nt,))
    @assert(size(contents["tmhadiv"]) == (nt,))
    return true
end

function test_populate_grid_ggd()
    filename = "samples/b2fgmtry"
    gmtry = read_b2_output(filename)
    nxg, nyg = gmtry["nx,ny"]
    nx = nxg + 2  # Account for guard cells
    ny = nyg + 2  # Account for guard cells
    filename = "samples/b2time.nc"
    b2t = read_b2_output(filename)
    crx = reshape(gmtry["crx"], (4, ny, nx))
    cry = reshape(gmtry["cry"], (4, ny, nx))
    populate_grid_ggd("samples/gridspacedesc.yml", crx, cry, "electrons", "density", b2t["ne2d"], b2t["timesa"])
    return true
end

@testset "omasstuff" begin
    @test try_omas() === nothing
    @test test_read_b2_output()
    @test test_populate_grid_ggd()
end
