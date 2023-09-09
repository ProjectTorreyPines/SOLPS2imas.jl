import SOLPS2IMAS
using Test
using YAML: load_file as YAML_load_file

function test_generate_test_data()
    cr, cz = generate_test_data()
    plot(cr, cz)
    return true
end

function test_ind_conversion()
    nx = 92
    ny = 38
    success = true
    for iy = 1:ny
        for ix = 1:nx
            cxy = SOLPS2IMAS.ctoxy(SOLPS2IMAS.xytoc(ix, iy; nx=nx); nx=nx)
            if cxy != (ix, iy)
                ic = SOLPS2IMAS.xytoc(ix, iy; nx=nx)
                println("ic: ", ic, ", (ix, iy): (", ix, ", ", iy, "), converted_ix,iy: ", cxy)
                success = false
            end
        end
    end

    for ic = 1:nx*ny
        cic = SOLPS2IMAS.xytoc(SOLPS2IMAS.ctoxy(ic; nx=nx)...; nx=nx)
        if cic != ic
            ix, iy = SOLPS2IMAS.ctoxy(ic; nx=nx)
            println("ic: ", ic, ", (ix, iy): (", ix, ", ", iy, "), converted_ic: ", cic)
            success = false
        end
    end
    return success
end

function test_read_b2_output()
    contents = SOLPS2IMAS.read_b2_output("$(@__DIR__)/../samples/b2fstate")
    nt = contents["dim"]["time"]
    nx = contents["dim"]["nx"]
    ny = contents["dim"]["ny"]
    ns = contents["dim"]["ns"]
    @assert(size(contents["data"]["te"]) == (nt, ny, nx))
    @assert(size(contents["data"]["na"]) == (nt, ns, ny, nx))
    @assert(size(contents["data"]["fna"]) == (nt, ns, 2, ny, nx))
    @assert(size(contents["data"]["fhe"]) == (nt, ns, ny, nx))

    contents = SOLPS2IMAS.read_b2_output("$(@__DIR__)/../samples/b2fgmtry")
    nt = contents["dim"]["time"]
    nxg = contents["dim"]["nx"]
    nyg = contents["dim"]["ny"]
    @assert(size(contents["data"]["crx"]) == (nt, 4, nyg, nxg))
    @assert(size(contents["data"]["cry"]) == (nt, 4, nyg, nxg))
    @assert(nyg == ny)
    @assert(nxg == nx)

    contents = SOLPS2IMAS.read_b2_output("$(@__DIR__)/../samples/b2time_red.nc")
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
    b2gmtry = "$(@__DIR__)/../samples/b2fgmtry"
    b2output = "$(@__DIR__)/../samples/b2time_red.nc"
    gsdesc = "$(@__DIR__)/../samples/gridspacedesc.yml"
    b2t = SOLPS2IMAS.read_b2_output(b2output)
    nx = b2t["dim"]["nx"]
    @time dd = SOLPS2IMAS.solps2imas(b2gmtry, b2output, gsdesc)
    # Check time stamp 3 at iy=4, ix=5
    it = 3
    iy = 4
    ix = 5
    @assert(b2t["data"]["ne2d"][3, iy, ix] == dd.edge_profiles.ggd[3].electrons.density[5].values[(iy - 1) * nx + ix])
    return true
end


@testset "omasstuff" begin
    @test SOLPS2IMAS.try_omas() === nothing
    @test test_ind_conversion()
    @test test_read_b2_output()
    @test test_solps2imas()
end

@testset "utilities" begin
    gsdesc = YAML_load_file("$(@__DIR__)/../samples/gridspacedesc.yml")
    println("hello there")
    @test SOLPS2IMAS.find_subset_index(gsdesc, 1) == 1
    @test SOLPS2IMAS.find_subset_index(gsdesc, 5) == 5
    @test SOLPS2IMAS.find_subset_index(gsdesc, 101) == 27
end
