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
    contents = read_b2_output("samples/b2fstate")
    nt = contents["dim"]["time"]
    nx = contents["dim"]["nx"]
    ny = contents["dim"]["ny"]
    ns = contents["dim"]["ns"]
    @assert(size(contents["data"]["te"]) == (nt, ny, nx))
    @assert(size(contents["data"]["na"]) == (nt, ns, ny, nx))
    @assert(size(contents["data"]["fna"]) == (nt, ns, 2, ny, nx))
    @assert(size(contents["data"]["fhe"]) == (nt, ns, ny, nx))

    contents = read_b2_output("samples/b2fgmtry")
    nt = contents["dim"]["time"]
    nxg = contents["dim"]["nx"]
    nyg = contents["dim"]["ny"]
    @assert(size(contents["data"]["crx"]) == (nt, 4, nyg, nxg))
    @assert(size(contents["data"]["cry"]) == (nt, 4, nyg, nxg))
    @assert(nyg == ny)
    @assert(nxg == nx)

    contents = read_b2_output("samples/b2time.nc")
    nx = contents["dim"]["nx"]
    ny = contents["dim"]["ny"]
    nybl = contents["dim"]["nybl"]
    nybr = contents["dim"]["nybr"]
    nt = contents["dim"]["time"]
    ns = contents["dim"]["ns"]
    @assert(size(contents["data"]["ne2d"]) == (nt, ny, nx))
    @assert(size(contents["data"]["ft3dl"]) == (nt, nybl))
    @assert(size(contents["data"]["fl3dl"]) == (nt, nybl))
    @assert(size(contents["data"]["fc3dl"]) == (nt, nybl))
    @assert(size(contents["data"]["ft3dr"]) == (nt, nybr))
    @assert(size(contents["data"]["fl3dr"]) == (nt, nybr))
    @assert(size(contents["data"]["fc3dr"]) == (nt, nybr))
    @assert(size(contents["data"]["fna3da"]) == (nt, ns, nybr))
    @assert(size(contents["data"]["tmhacore"]) == (nt, ))
    @assert(size(contents["data"]["tmhasol"]) == (nt,))
    @assert(size(contents["data"]["tmhadiv"]) == (nt,))
    return true
end

function test_solps2imas()
    b2gmtry = "samples/b2fgmtry"
    b2output = "samples/b2time.nc"
    gsdesc = "samples/gridspacedesc.yml"
    b2t = read_b2_output(b2output)
    nx = b2t["dim"]["nx"]
    dd = solps2imas(b2gmtry, b2output, gsdesc)
    # Check time stamp 3 at iy=4, ix=5
    it = 3
    iy = 4
    ix = 5
    @assert(b2t["data"]["ne2d"][3, iy, ix] == dd.edge_profiles.ggd[3].electrons.density[5].values[(iy - 1) * nx + ix])
    return true
end

function test_populate_grid_ggd()
    gmtry = read_b2_output("samples/b2fgmtry")
    nx = gmtry["dim"]["nx"]
    ny = gmtry["dim"]["ny"]
    b2t = read_b2_output("samples/b2time.nc")
    crx = reshape(gmtry["data"]["crx"], (4, ny, nx))
    cry = reshape(gmtry["data"]["cry"], (4, ny, nx))
    populate_grid_ggd("samples/gridspacedesc.yml", crx, cry, "electrons", "density", b2t["data"]["ne2d"], b2t["data"]["timesa"])
    return true
end

@testset "omasstuff" begin
    @test try_omas() === nothing
    @test test_read_b2_output()
    @test test_solps2imas()
end
