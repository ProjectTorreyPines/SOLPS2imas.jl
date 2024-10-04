using SOLPS2imas: SOLPS2imas
using Test
using YAML: load_file as YAML_load_file
using ArgParse: ArgParse
using IMASdd: IMASdd
import SOLPS2imas: get_grid_subset # , read_b2_boundary_parameters

allowed_rtol = 1e-4

# To run a subset of tests while executing via > include("test/runtests.jl"), define newARGS. For example:
#   newARGS = String["--ind", "--b2]
# to run just the ind and b2 test sets.

function parse_commandline()
    localARGS = @isdefined(newARGS) ? newARGS : ARGS  # Thanks https://stackoverflow.com/a/44978474/6605826
    s = ArgParse.ArgParseSettings(; description="Run tests. Default is all tests.")

    ArgParse.add_arg_table!(s,
        ["--ind"],
        Dict(:help => "Test index conversions",
            :action => :store_true),
        ["--b2"],
        Dict(:help => "Test read_b2_output()",
            :action => :store_true),
        ["--solps2imas"],
        Dict(:help => "Test solps2imas (overall workflow)",
            :action => :store_true),
        ["--parser"],
        Dict(:help => "Stress test file parsing (other than b2 output files)",
            :action => :store_true),
        ["--fort"],
        Dict(:help => "Test triangular mesh generation from fort files",
            :action => :store_true),
        # ["--namelist"],
        # Dict(:help => "Test parsing of namelists",
        #     :action => :store_true),
    )
    args = ArgParse.parse_args(localARGS, s)
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
                cxy = SOLPS2imas.ctoxy(SOLPS2imas.xytoc(ix, iy; nx=nx); nx=nx)
                if cxy != (ix, iy)
                    ic = SOLPS2imas.xytoc(ix, iy; nx=nx)
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
            cic = SOLPS2imas.xytoc(SOLPS2imas.ctoxy(ic; nx=nx)...; nx=nx)
            if cic != ic
                ix, iy = SOLPS2imas.ctoxy(ic; nx=nx)
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

if args["parser"]
    @testset "Test file parsing in depth" begin
        b2mn_samples = "$(@__DIR__)/../samples/" .* [
            "b2mn.dat",
            "test_b2mn.dat",
        ]
        always_required_keys = ["b2mndr_ntim", "b2mndr_dtim"]
        for b2mn_sample ∈ b2mn_samples
            b2mn_data = SOLPS2imas.read_b2mn_output(b2mn_sample)
            for ark ∈ always_required_keys
                @test ark in keys(b2mn_data)
            end
            b2mn_json = IMASdd.JSON.parsefile(b2mn_sample * ".json")
            @test b2mn_json == b2mn_data
        end
    end
end

if args["b2"]
    @testset "Test read_b2_output" begin
        contents = SOLPS2imas.read_b2_output("$(@__DIR__)/../samples/b2fstate")
        nt = contents["dim"]["time"]
        nx = contents["dim"]["nx"]
        ny = contents["dim"]["ny"]
        ns = contents["dim"]["ns"]
        @test size(contents["data"]["te"]) == (nt, ny, nx)
        @test size(contents["data"]["na"]) == (nt, ns, ny, nx)
        @test size(contents["data"]["fna"]) == (nt, ns, 2, ny, nx)
        @test size(contents["data"]["fhe"]) == (nt, ns, ny, nx)

        contents = SOLPS2imas.read_b2_output("$(@__DIR__)/../samples/b2fgmtry")
        nt = contents["dim"]["time"]
        nxg = contents["dim"]["nx"]
        nyg = contents["dim"]["ny"]
        @test size(contents["data"]["crx"]) == (nt, 4, nyg, nxg)
        @test size(contents["data"]["cry"]) == (nt, 4, nyg, nxg)
        @test nyg == ny
        @test nxg == nx

        contents = SOLPS2imas.read_b2_output("$(@__DIR__)/../samples/b2time_red.nc")
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
    @testset "Test solps2imas() (overall workflow)" begin
        b2gmtry = "$(@__DIR__)/../samples/b2fgmtry"
        b2output = "$(@__DIR__)/../samples/b2time.nc"
        b2mn = "$(@__DIR__)/../samples/b2mn.dat"
        b2t = SOLPS2imas.read_b2_output(b2output)
        nx = b2t["dim"]["nx"]
        print("solps2imas() time: ")
        @time dd = SOLPS2imas.solps2imas(b2gmtry, b2output; b2mn=b2mn, load_bb=false)
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
        gmtry = SOLPS2imas.read_b2_output(b2gmtry)
        ny = gmtry["dim"]["ny"]
        cuts = Dict([(Symbol(key), gmtry["data"][key][1]) for key ∈ cut_keys])
        subset_pfrcut =
            SOLPS2imas.get_grid_subset(dd.edge_profiles.grid_ggd[1], 8)
        subset_corebnd =
            SOLPS2imas.get_grid_subset(dd.edge_profiles.grid_ggd[1], 15)
        subset_separatrix =
            SOLPS2imas.get_grid_subset(dd.edge_profiles.grid_ggd[1], 16)
        cells = dd.edge_profiles.grid_ggd[1].space[1].objects_per_dimension[3].object
        subset_pfrcut_element_list =
            [ele.object[1].index for ele ∈ subset_pfrcut.element]
        subset_corebnd_element_list =
            [ele.object[1].index for ele ∈ subset_corebnd.element]
        subset_separatrix_element_list =
            [ele.object[1].index for ele ∈ subset_separatrix.element]
        brute_force_pfrcut_list = []
        brute_force_corebnd_list = []
        brute_force_separatrix_list = []
        for iy ∈ 1:ny
            for ix ∈ 1:nx
                for boundary_ind ∈ 1:4
                    edge_ind =
                        cells[SOLPS2imas.xytoc(ix, iy; nx)].boundary[boundary_ind].index
                    if SOLPS2imas.is_pfr_cut(; ix, iy, cells, nx, boundary_ind, cuts...)
                        append!(brute_force_pfrcut_list, edge_ind)
                    elseif SOLPS2imas.is_core_boundary(; ix, iy, boundary_ind, cuts...)
                        append!(brute_force_corebnd_list, edge_ind)
                    elseif SOLPS2imas.is_separatrix(; iy, boundary_ind, cuts...)
                        append!(brute_force_separatrix_list, edge_ind)
                    end
                end
            end
        end
        @test Set(brute_force_pfrcut_list) == Set(subset_pfrcut_element_list)
        @test Set(brute_force_corebnd_list) == Set(subset_corebnd_element_list)
        @test Set(brute_force_separatrix_list) == Set(subset_separatrix_element_list)

        # Test loading of magnetic field data from b2fgmtry to equilibrium IDS
        @time ddbb = SOLPS2imas.solps2imas(b2gmtry, b2output; b2mn=b2mn, load_bb=true)
        btor = ddbb.equilibrium.time_slice[1].ggd[1].b_field_tor[1].values
        bz = ddbb.equilibrium.time_slice[1].ggd[1].b_field_z[1].values
        br = ddbb.equilibrium.time_slice[1].ggd[1].b_field_r[1].values
        @test length(btor) == (nx * ny)
        @test length(bz) == length(btor)
        @test length(br) == length(bz)
    end
end

if args["fort"]
    @testset "Test triangular mesh generation from fort files" begin
        fort = (
            "$(@__DIR__)/../samples/fort.33",
            "$(@__DIR__)/../samples/fort.34",
            "$(@__DIR__)/../samples/fort.35")
        b2gmtry = "$(@__DIR__)/../samples/b2fgmtry"
        b2output = "$(@__DIR__)/../samples/b2time.nc"
        b2mn = "$(@__DIR__)/../samples/b2mn.dat"
        ids = SOLPS2imas.solps2imas(b2gmtry, b2output; b2mn=b2mn, fort=fort)
        grid_ggd = ids.edge_profiles.grid_ggd[1]
        space = grid_ggd.space[1]

        subset_nodes = get_grid_subset(grid_ggd, 1)
        subset_faces = get_grid_subset(grid_ggd, "faces")
        subset_cells = get_grid_subset(grid_ggd, "cells")
        subset_b25nodes = get_grid_subset(grid_ggd, -1)
        subset_b25faces = get_grid_subset(grid_ggd, -2)
        subset_b25cells = get_grid_subset(grid_ggd, -5)
        subset_trinodes = get_grid_subset(grid_ggd, -101)
        subset_trifaces = get_grid_subset(grid_ggd, -102)
        subset_tricells = get_grid_subset(grid_ggd, -105)
        subset_comnodes = get_grid_subset(grid_ggd, -201)
        subset_comfaces = get_grid_subset(grid_ggd, -202)
        subset_comcells = get_grid_subset(grid_ggd, -205)

        nodes = grid_ggd.space[1].objects_per_dimension[1].object
        edges = grid_ggd.space[1].objects_per_dimension[2].object
        cells = grid_ggd.space[1].objects_per_dimension[3].object

        gmtry = SOLPS2imas.read_b2_output(b2gmtry)
        nx = gmtry["dim"]["nx"]
        ny = gmtry["dim"]["ny"]
        extra_nodes = 2 * (nx + ny)

        @test length(subset_nodes.element) == length(nodes)
        @test length(subset_faces.element) == length(edges)
        @test length(subset_cells.element) == length(cells)

        @test length(subset_b25nodes.element) ==
              length(subset_comnodes.element) + extra_nodes

        @test length(subset_trinodes.element) + extra_nodes == length(nodes)

        for ele ∈ subset_tricells.element
            tricell_ind = ele.object[1].index
            for bnd_ind ∈ 1:3
                if length(cells[tricell_ind].boundary[bnd_ind].neighbours) > 0
                    common_nodes = intersect(
                        cells[tricell_ind].nodes,
                        cells[cells[tricell_ind].boundary[bnd_ind].neighbours[1]].nodes,
                    )
                    vert_no = Tuple(indexin(common_nodes, cells[tricell_ind].nodes))
                    @test SOLPS2imas.chosen_tri_edge_order[bnd_ind][2] == vert_no
                end
            end
        end
    end
end

# if args["namelist"]
#     @testset "Test parsing of namelists" begin
#         # Basic parameters namelist parsing
#         testfile = "$(@__DIR__)/../samples/b2.boundary.parameters"
#         boundary_params = SOLPS2imas.read_b2_boundary_parameters(testfile)
#         println(boundary_params)
#         @test boundary_params["power_electrons"] > 0.0
#         @test boundary_params["power_ions"] > 0.0
#         @test boundary_params["number_of_boundaries"] >
#               boundary_params["number_of_core_source_boundaries"]

#         # Using parameters namelist to populate summary data
#         ids = IMASdd.dd()
#         SOLPS2imas.load_summary_data!(ids, (testfile, "", "", ""))
#         @test !(ismissing(ids.summary.heating_current_drive.power_ec, :value))
#     end
# end
