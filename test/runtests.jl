using SOLPS2IMAS: SOLPS2IMAS
using Test
using YAML: load_file as YAML_load_file
using ArgParse: ArgParse

allowed_rtol = 1e-4

function parse_commandline()
    s = ArgParse.ArgParseSettings(; description="Run tests. Default is all tests.")

    ArgParse.add_arg_table!(s,
        ["--ind"],
        Dict(:help => "Test index conversions",
            :action => :store_true),
        ["--b2"],
        Dict(:help => "Test read_b2_output()",
            :action => :store_true),
        ["--solps2imas"],
        Dict(:help => "Test solps2imas",
            :action => :store_true),
    )
    args = ArgParse.parse_args(s)
    if !any(values(args)) # If no flags are set, run all tests
        for k ∈ keys(args)
            args[k] = true
        end
    end
    return args
end
args = parse_commandline()

if args["ind"]
    @testset "Test index conversions" begin
        nx = 92
        ny = 38
        success = true
        for iy ∈ 1:ny
            for ix ∈ 1:nx
                cxy = SOLPS2IMAS.ctoxy(SOLPS2IMAS.xytoc(ix, iy; nx=nx); nx=nx)
                if cxy != (ix, iy)
                    ic = SOLPS2IMAS.xytoc(ix, iy; nx=nx)
                    println(
                        "ic: ",
                        ic,
                        ", (ix, iy): (",
                        ix,
                        ", ",
                        iy,
                        "), converted_ix,iy: ",
                        cxy,
                    )
                    success = false
                end
            end
        end

        for ic ∈ 1:nx*ny
            cic = SOLPS2IMAS.xytoc(SOLPS2IMAS.ctoxy(ic; nx=nx)...; nx=nx)
            if cic != ic
                ix, iy = SOLPS2IMAS.ctoxy(ic; nx=nx)
                println(
                    "ic: ",
                    ic,
                    ", (ix, iy): (",
                    ix,
                    ", ",
                    iy,
                    "), converted_ic: ",
                    cic,
                )
                success = false
            end
        end
        @test success
    end
end

if args["b2"]
    @testset "Test read_b2_output" begin
        contents = SOLPS2IMAS.read_b2_output("$(@__DIR__)/../samples/b2fstate")
        nt = contents["dim"]["time"]
        nx = contents["dim"]["nx"]
        ny = contents["dim"]["ny"]
        ns = contents["dim"]["ns"]
        @test size(contents["data"]["te"]) == (nt, ny, nx)
        @test size(contents["data"]["na"]) == (nt, ns, ny, nx)
        @test size(contents["data"]["fna"]) == (nt, ns, 2, ny, nx)
        @test size(contents["data"]["fhe"]) == (nt, ns, ny, nx)

        contents = SOLPS2IMAS.read_b2_output("$(@__DIR__)/../samples/b2fgmtry")
        nt = contents["dim"]["time"]
        nxg = contents["dim"]["nx"]
        nyg = contents["dim"]["ny"]
        @test size(contents["data"]["crx"]) == (nt, 4, nyg, nxg)
        @test size(contents["data"]["cry"]) == (nt, 4, nyg, nxg)
        @test nyg == ny
        @test nxg == nx

        contents = SOLPS2IMAS.read_b2_output("$(@__DIR__)/../samples/b2time_red.nc")
        nx = contents["dim"]["nx"]
        ny = contents["dim"]["ny"]
        nybl = contents["dim"]["nybl"]
        nybr = contents["dim"]["nybr"]
        nt = contents["dim"]["time"]
        ns = contents["dim"]["ns"]
        @test size(contents["data"]["ne2d"]) == (nt, ny, nx)
        @test size(contents["data"]["ft3dl"]) == (nt, nybl)
        @test size(contents["data"]["fl3dl"]) == (nt, nybl)
        @test size(contents["data"]["fc3dl"]) == (nt, nybl)
        @test size(contents["data"]["ft3dr"]) == (nt, nybr)
        @test size(contents["data"]["fl3dr"]) == (nt, nybr)
        @test size(contents["data"]["fc3dr"]) == (nt, nybr)
        @test size(contents["data"]["fna3da"]) == (nt, ns, nybr)
        @test size(contents["data"]["tmhacore"]) == (nt,)
        @test size(contents["data"]["tmhasol"]) == (nt,)
        @test size(contents["data"]["tmhadiv"]) == (nt,)
    end
end

if args["solps2imas"]
    @testset "Test solps2imas()" begin
        b2gmtry = "$(@__DIR__)/../samples/b2fgmtry"
        b2output = "$(@__DIR__)/../samples/b2time.nc"
        gsdesc = "$(@__DIR__)/../samples/gridspacedesc.yml"
        b2mn = "$(@__DIR__)/../samples/b2mn.dat"
        b2t = SOLPS2IMAS.read_b2_output(b2output)
        nx = b2t["dim"]["nx"]
        print("solps2imas() time: ")
        @time dd = SOLPS2IMAS.solps2imas(b2gmtry, b2output, gsdesc, b2mn)
        # Check time stamp 3 at iy=4, ix=5
        it = 3
        iy = 4
        ix = 5
        @test (
            b2t["data"]["ne2d"][3, iy, ix]
            ==
            dd.edge_profiles.ggd[3].electrons.density[1].values[(iy-1)*nx+ix]
        )

        # Checking if subsets obtained from set operations are same
        # as using a brute force definition which is too dependent
        # on the correct ordering of nodes in SOLPS data files.
        cut_keys = ["leftcut", "rightcut", "bottomcut", "topcut"]
        gmtry = SOLPS2IMAS.read_b2_output(b2gmtry)
        ny = gmtry["dim"]["ny"]
        cuts = Dict([(Symbol(key), gmtry["data"][key][1]) for key ∈ cut_keys])
        subset_pfrcut =
            SOLPS2IMAS.get_grid_subset_with_index(dd.edge_profiles.grid_ggd[1], 8)
        subset_corebnd =
            SOLPS2IMAS.get_grid_subset_with_index(dd.edge_profiles.grid_ggd[1], 15)
        subset_separatix =
            SOLPS2IMAS.get_grid_subset_with_index(dd.edge_profiles.grid_ggd[1], 16)
        cells = dd.edge_profiles.grid_ggd[1].space[1].objects_per_dimension[3].object
        subset_pfrcut_element_list =
            [ele.object[1].index for ele ∈ subset_pfrcut.element]
        subset_corebnd_element_list =
            [ele.object[1].index for ele ∈ subset_corebnd.element]
        subset_separatix_element_list =
            [ele.object[1].index for ele ∈ subset_separatix.element]
        brute_force_pfrcut_list = []
        brute_force_corebnd_list = []
        brute_force_separatix_list = []
        for iy ∈ 1:ny
            for ix ∈ 1:nx
                for boundary_ind ∈ 1:4
                    edge_ind =
                        cells[SOLPS2IMAS.xytoc(ix, iy; nx)].boundary[boundary_ind].index
                    if SOLPS2IMAS.is_pfr_cut(; ix, iy, cells, nx, boundary_ind, cuts...)
                        append!(brute_force_pfrcut_list, edge_ind)
                    elseif SOLPS2IMAS.is_core_boundary(; ix, iy, boundary_ind, cuts...)
                        append!(brute_force_corebnd_list, edge_ind)
                    elseif SOLPS2IMAS.is_separatix(; iy, boundary_ind, cuts...)
                        append!(brute_force_separatix_list, edge_ind)
                    end
                end
            end
        end
        @test Set(brute_force_pfrcut_list) == Set(subset_pfrcut_element_list)
        @test Set(brute_force_corebnd_list) == Set(subset_corebnd_element_list)
        @test Set(brute_force_separatix_list) == Set(subset_separatix_element_list)
    end
end
