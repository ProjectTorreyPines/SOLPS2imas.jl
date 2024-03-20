module SOLPS2IMAS

using Revise
import OMAS as IMASDD
using NCDatasets: Dataset, dimnames
using YAML: load_file as YAML_load_file
using DelimitedFiles: readdlm
import GGDUtils: add_subset_element!, get_grid_subset_with_index, get_subset_boundary,
    get_subset_space, subset_do

export read_b2_output
export solps2imas

include("parser.jl")

include("subset_id.jl")

include("solps_mesh_tools.jl")

default_gsdesc = "$(@__DIR__)/default_grid_space_description.yml"

"""
    dict2prop!(obj, dict)

Copies grid_ggd and space description in dict format to the data structure recursively.
"""
dict2prop!(obj, dict) =
    for (key, prop) ∈ dict
        if isa(key, Int)
            if length(obj) < key
                resize!(obj, key)
            end
            dict2prop!(obj[key], prop)
        elseif isa(key, String)
            if isa(prop, Dict)
                dict2prop!(getfield(obj, Symbol(key)), prop)
            else
                setproperty!(obj, Symbol(key), prop)
            end
        end
    end

solps_var_to_imas = YAML_load_file("$(@__DIR__)/solps_var_to_imas.yml")
"""
    val_obj(
    ggd::IMASDD.edge_profiles__ggd,
    var::String,
    grid_ggd_index::Int64;
    gsi_ch::Dict{Int64, Int64}=Dict{Int64, Int64}(),

)

Given SOLPS variable name (var), returns pair of parent object and property name
to write value on. If var is not found in solps_var_to_imas, returns nothing, nothing.
Optionally, a mapping of possibly changes in grid_subset_index can be provided as
a dictionary gsi_ch.
"""
function val_obj(
    ggd::IMASDD.edge_profiles__ggd,
    var::String,
    grid_ggd_index::Int64;
    gsi_ch::Dict{Int64, Int64}=Dict{Int64, Int64}(),
)
    if var ∉ keys(solps_var_to_imas)
        return nothing, nothing
    else
        path, gsi = solps_var_to_imas[var]
        if gsi ∈ gsi_ch.keys
            gsi = gsi_ch[gsi]
        end
        parent = ggd
        path_fields = split(path, ".")
        for pf ∈ path_fields[1:end-1]
            if occursin("[", pf)
                parent = getfield(parent, Symbol(pf[1:findfirst('[', pf)-1]))
                ind_str = pf[findfirst('[', pf)+1:findfirst(']', pf)-1]
                if ind_str == ":"
                    resize!(parent, length(parent) + 1)
                    parent = parent[end]
                else
                    ind = parse(Int64, ind_str)
                    if length(parent) < ind
                        resize!(parent, ind)
                    end
                    parent = parent[ind]
                end
            else
                parent = getfield(parent, Symbol(pf))
            end
            if :grid_subset_index ∈ fieldnames(typeof(parent))
                parent.grid_subset_index = gsi
                parent.grid_index = grid_ggd_index
            end
        end
        return parent, Symbol(path_fields[end])
    end
end

# Following convention is used to index the edges of a cell
# This ends up going around the cell starting with bottom x-edge,
# right y-edge, top x-edge, and left y-edge
# Thus, x-edges will have odd boundary index and y_edges will have even
# List of tuples (boundary_ind, (corner pair forming edge))
chosen_edge_order = [(1, (1, 2)),
    (2, (2, 4)),
    (3, (4, 3)),
    (4, (3, 1))]

# Following convention is used to index the edges of a triangular cell
chosen_tri_edge_order = [(1, (1, 2)),
    (2, (2, 3)),
    (3, (1, 3))]

"""
    solps2imas(
    b2gmtry::String,
    b2output::String;
    gsdesc::String=default_gsdesc,
    b2mn::String="",
    fort::Tuple{String, String, String}=("", "", ""),
    fort_tol::Float64=1e-6,
    load_bb::Bool=false,

)

Main function of the module. Takes in a geometry file and a
output file (either b2time or b2fstate) and a grid_ggd
description in the form of a Dict or filename to equivalent
YAML file. Additionally, EIRENE fort files can be provided as tuple of 3 filenames
consisting fort.33, fort.34, and fort.35 files. The grids in these files are matched
with SOLPS grid with a tolerance of fort_tol (defaults to 1e-6).
Returns data in IMASDD.dd datastructure.
"""
function solps2imas(
    b2gmtry::String,
    b2output::String;
    gsdesc::String=default_gsdesc,
    b2mn::String="",
    fort::Tuple{String, String, String}=("", "", ""),
    fort_tol::Float64=1e-6,
    load_bb::Bool=false,
)
    # Initialize an empty IMAS data structre
    ids = IMASDD.dd()

    # Setup the grid first
    gmtry = read_b2_output(b2gmtry)

    jxi = jxa = nothing
    if b2mn != ""
        mn = read_b2mn_output(b2mn)
        if "b2mwti_jxa" ∈ keys(mn)
            jxa = mn["b2mwti_jxa"]
        end
        if "b2mwti_jxi" ∈ keys(mn)
            jxi = mn["b2mwti_jxi"]
        end
    end

    nx = gmtry["dim"]["nx"]
    ny = gmtry["dim"]["ny"]
    ncell = nx * ny
    crx = gmtry["data"]["crx"]
    cry = gmtry["data"]["cry"]
    cut_keys = ["leftcut", "rightcut", "bottomcut", "topcut"]
    cuts_found = cut_keys ⊆ keys(gmtry["data"])
    if cuts_found
        cuts = Dict([(Symbol(key), gmtry["data"][key][1]) for key ∈ cut_keys])
        xpoints_nodes = get_xpoint_nodes(gmtry)
    end

    if typeof(gsdesc) == String
        gsdesc = YAML_load_file(gsdesc)
    end

    # Add grid_ggd array equal to number of time steps
    resize!(ids.edge_profiles.grid_ggd, gmtry["dim"]["time"])
    for it ∈ 1:gmtry["dim"]["time"]
        grid_ggd = ids.edge_profiles.grid_ggd[it]
        grid_ggd.time = Float64.(gmtry["data"]["timesa"][it])
        dict2prop!(grid_ggd, gsdesc)
        for sn ∈ keys(gsdesc["space"])
            space = grid_ggd.space[sn]
            # Assuming following to be standard for now.
            # We can add this info through YAML as well
            resize!(space.objects_per_dimension, 3)
            o1 = space.objects_per_dimension[1]  # 1D objects
            o2 = space.objects_per_dimension[2]  # 2D objects
            o3 = space.objects_per_dimension[3]  # 3D objects

            subset_nodes = get_grid_subset_with_index(grid_ggd, 1)
            subset_faces = get_grid_subset_with_index(grid_ggd, 2)
            subset_xedges = get_grid_subset_with_index(grid_ggd, 3)
            subset_yedges = get_grid_subset_with_index(grid_ggd, 4)
            subset_cells = get_grid_subset_with_index(grid_ggd, 5)
            if cuts_found
                subset_core = get_grid_subset_with_index(grid_ggd, 22)
                subset_sol = get_grid_subset_with_index(grid_ggd, 23)
                subset_odr = get_grid_subset_with_index(grid_ggd, 24)
                subset_idr = get_grid_subset_with_index(grid_ggd, 25)
                subset_xp = get_grid_subset_with_index(grid_ggd, 6)
                subset_corecut = get_grid_subset_with_index(grid_ggd, 7)
                subset_pfrcut = get_grid_subset_with_index(grid_ggd, 8)
                subset_othroat = get_grid_subset_with_index(grid_ggd, 9)
                subset_ithroat = get_grid_subset_with_index(grid_ggd, 10)
                if !isnothing(jxa)
                    subset_omp = get_grid_subset_with_index(grid_ggd, 11)
                    subset_ompsep = get_grid_subset_with_index(grid_ggd, 101)
                end
                if !isnothing(jxi)
                    subset_imp = get_grid_subset_with_index(grid_ggd, 12)
                    subset_impsep = get_grid_subset_with_index(grid_ggd, 102)
                end
                subset_corebnd = get_grid_subset_with_index(grid_ggd, 15)
                subset_separatix = get_grid_subset_with_index(grid_ggd, 16)
                subset_otsep = get_grid_subset_with_index(grid_ggd, 103)
                subset_itsep = get_grid_subset_with_index(grid_ggd, 104)
            end

            subset_otarget = get_grid_subset_with_index(grid_ggd, 13)
            subset_itarget = get_grid_subset_with_index(grid_ggd, 14)

            # Resizing objects to hold cell geometry data
            # Should be fewer than this many points, but this way we won't under-fill
            # nodes = resize!(o1.object, ncell * 4)  # Nodes (1D)
            # edges = resize!(o2.object, ncell * 4)  # Edges (2D)
            nodes = o1.object                  # Nodes (1D)
            edges = o2.object                  # Edges (2D)
            cells = resize!(o3.object, ncell)  # Cells (3D)

            # Initialize geometry for 1D objects(nodes), nodes for 2D objects(edges)
            # for i ∈ 1:(ncell*4)
            #     nodes[i].geometry = [0.0, 0.0]
            #     edges[i].nodes = [0, 0]
            #     resize!(edges[i].boundary, 2)
            #     for bnd ∈ edges[i].boundary
            #         bnd.neighbours = Int64[]
            #     end
            # end
            # Initialize nodes and boundaries for cells
            for i ∈ 1:(ncell)
                cells[i].nodes = [0, 0, 0, 0]
                resize!(cells[i].boundary, 4)
                for bnd ∈ cells[i].boundary
                    bnd.neighbours = Int64[]
                end
            end

            j = 1
            edge_ind = 1
            # Setting up space with nodes, edges and cells
            for iy ∈ 1:ny
                for ix ∈ 1:nx
                    ic::Int = (iy - 1) * nx + ix
                    # Adding node positions and cell corners data
                    for icorner ∈ 1:4
                        # Have to search to see if the node is already added and then
                        # record its index
                        # If not already listed, then list it under new index and
                        # record that
                        # Note that time index has been fixed to 1 here. Only handling
                        # fixed grid geometry through the run cases.
                        i_existing = search_points(
                            nodes,
                            crx[1, icorner, iy, ix],
                            cry[1, icorner, iy, ix],
                        )[1]
                        if i_existing == 0
                            resize!(nodes, j)
                            nodes[j].geometry =
                                [crx[1, icorner, iy, ix], cry[1, icorner, iy, ix]]
                            cells[ic].nodes[icorner] = j
                            add_subset_element!(subset_nodes, sn, 1, j)
                            if cuts_found && xpoints_nodes[it][1] == nodes[j].geometry
                                add_subset_element!(subset_xp, sn, 1, j)
                            end
                            j += 1
                        else
                            cells[ic].nodes[icorner] = i_existing[1]
                        end
                    end
                    # Adding edges (faces with toroidal elongation)
                    # Adding same edges as boundary to cells
                    for (boundary_ind, edge_pair) ∈ chosen_edge_order
                        edge_nodes = [cells[ic].nodes[icorner] for icorner ∈ edge_pair]
                        existing_edge_ind = search_edges(edges, edge_nodes)
                        if existing_edge_ind == 0
                            resize!(edges, edge_ind)
                            edges[edge_ind].nodes = edge_nodes
                            resize!(edges[edge_ind].boundary, 2)
                            for (ii, edge_bnd) ∈ enumerate(edges[edge_ind].boundary)
                                edge_bnd.index = edge_nodes[ii]
                            end
                            edges[edge_ind].measure =
                                distance_between_nodes(nodes, edge_nodes)
                            cells[ic].boundary[boundary_ind].index = edge_ind
                            add_subset_element!(subset_faces, sn, 2, edge_ind)
                            add_subset_element!(
                                subset_xedges,
                                sn,
                                2,
                                edge_ind,
                                is_x_aligned;
                                boundary_ind,
                            )
                            add_subset_element!(
                                subset_yedges,
                                sn,
                                2,
                                edge_ind,
                                is_y_aligned;
                                boundary_ind,
                            )
                            edge_ind += 1
                        else
                            cells[ic].boundary[boundary_ind].index = existing_edge_ind
                        end
                    end
                    add_subset_element!(subset_cells, sn, 3, ic)
                    if cuts_found
                        add_subset_element!(
                            subset_core,
                            sn,
                            3,
                            ic,
                            in_core;
                            ix,
                            iy,
                            cuts...,
                        )
                        add_subset_element!(subset_sol, sn, 3, ic, in_sol; iy, cuts...)
                        add_subset_element!(
                            subset_idr,
                            sn,
                            3,
                            ic,
                            in_idr;
                            ix,
                            iy,
                            cuts...,
                        )
                        add_subset_element!(
                            subset_odr,
                            sn,
                            3,
                            ic,
                            in_odr;
                            ix,
                            iy,
                            cuts...,
                        )
                    end
                end
            end
            if cuts_found
                # Add boundaries
                attach_neightbours(cells, edges, gmtry, it)
                # Adding edges to subsets
                for iy ∈ 1:ny
                    for ix ∈ 1:nx
                        for boundary_ind ∈ 1:4
                            edge_ind =
                                cells[xytoc(ix, iy; nx)].boundary[boundary_ind].index
                            add_subset_element!(
                                subset_corecut,
                                sn,
                                2,
                                edge_ind,
                                is_core_cut;
                                ix,
                                iy,
                                cells,
                                nx,
                                boundary_ind,
                                cuts...,
                            )
                            if !isnothing(jxa)
                                add_subset_element!(
                                    subset_omp,
                                    sn,
                                    2,
                                    edge_ind,
                                    is_outer_midplane;
                                    ix,
                                    iy,
                                    jxa,
                                    boundary_ind,
                                    cuts...,
                                )
                            end
                            if !isnothing(jxi)
                                add_subset_element!(
                                    subset_imp,
                                    sn,
                                    2,
                                    edge_ind,
                                    is_inner_midplane;
                                    ix,
                                    iy,
                                    jxi,
                                    boundary_ind,
                                    cuts...,
                                )
                            end
                            add_subset_element!(
                                subset_othroat,
                                sn,
                                2,
                                edge_ind,
                                is_outer_throat;
                                ix,
                                iy,
                                boundary_ind,
                                cuts...,
                            )
                            add_subset_element!(
                                subset_ithroat,
                                sn,
                                2,
                                edge_ind,
                                is_inner_throat;
                                ix,
                                iy,
                                boundary_ind,
                                cuts...,
                            )
                            add_subset_element!(
                                subset_otarget,
                                sn,
                                2,
                                edge_ind,
                                is_outer_target;
                                ix,
                                nx,
                                boundary_ind,
                            )
                            add_subset_element!(
                                subset_itarget,
                                sn,
                                2,
                                edge_ind,
                                is_inner_target;
                                ix,
                                boundary_ind,
                            )
                        end
                    end
                end
                core_boundary_elements = get_subset_boundary(space, subset_core)
                sol_boundary_elements = get_subset_boundary(space, subset_sol)
                idr_boundary_elements = get_subset_boundary(space, subset_idr)
                odr_boundary_elements = get_subset_boundary(space, subset_odr)
                subset_pfrcut.element =
                    subset_do(intersect, idr_boundary_elements, odr_boundary_elements)
                subset_corebnd.element =
                    subset_do(setdiff, core_boundary_elements, sol_boundary_elements)
                subset_separatix.element = subset_do(intersect, sol_boundary_elements,
                    subset_do(union, core_boundary_elements,
                        odr_boundary_elements,
                        idr_boundary_elements))
                if !isnothing(jxa)
                    subset_ompsep.element = subset_do(
                        intersect,
                        subset_separatix.element,
                        subset_omp.element;
                        space,
                        use_nodes=true,
                    )
                end
                if !isnothing(jxi)
                    subset_impsep.element = subset_do(
                        intersect,
                        subset_separatix.element,
                        subset_imp.element;
                        space,
                        use_nodes=true,
                    )
                end
                subset_otsep.element = subset_do(
                    intersect,
                    subset_separatix.element,
                    subset_otarget.element;
                    space,
                    use_nodes=true,
                )
                subset_itsep.element = subset_do(
                    intersect,
                    subset_separatix.element,
                    subset_itarget.element;
                    space,
                    use_nodes=true,
                )
            end
        end  # End of setting up space
    end

    gsi_ch = Dict{Int64, Int64}()
    if fort != ("", "", "")
        grid_ggd = ids.edge_profiles.grid_ggd[1]
        space = grid_ggd.space[1]
        o1 = space.objects_per_dimension[1]  # 1D objects
        o2 = space.objects_per_dimension[2]  # 2D objects
        o3 = space.objects_per_dimension[3]  # 3D objects
        nodes = o1.object  # Nodes (1D)
        edges = o2.object  # Edges (2D)
        cells = o3.object  # Cells (3D)
        cur_no_subsets = length(grid_ggd.grid_subset)
        # Add 9 more subsets
        resize!(grid_ggd.grid_subset, cur_no_subsets + 12)
        # Copy all B2.5 nodes, faces, edges, cells to grid_subset with negative indices
        for (ii, gsi) ∈ enumerate([1, 2, 5])
            grid_ggd.grid_subset[cur_no_subsets+ii] =
                deepcopy(get_grid_subset_with_index(grid_ggd, gsi))
            grid_ggd.grid_subset[cur_no_subsets+ii].identifier.index = -gsi
            grid_ggd.grid_subset[cur_no_subsets+ii].identifier.name *= "_B2.5"
            gsi_ch[gsi] = -gsi
        end
        subset_nodes = get_grid_subset_with_index(grid_ggd, 1)
        subset_faces = get_grid_subset_with_index(grid_ggd, 2)
        subset_cells = get_grid_subset_with_index(grid_ggd, 5)
        subset_b25nodes = grid_ggd.grid_subset[cur_no_subsets+1]
        subset_b25faces = grid_ggd.grid_subset[cur_no_subsets+2]
        subset_b25cells = grid_ggd.grid_subset[cur_no_subsets+3]
        subset_trinodes = grid_ggd.grid_subset[cur_no_subsets+4]
        subset_trifaces = grid_ggd.grid_subset[cur_no_subsets+5]
        subset_tricells = grid_ggd.grid_subset[cur_no_subsets+6]
        subset_comnodes = grid_ggd.grid_subset[cur_no_subsets+7]
        subset_comfaces = grid_ggd.grid_subset[cur_no_subsets+8]
        subset_comcells = grid_ggd.grid_subset[cur_no_subsets+9]
        subset_extnodes = grid_ggd.grid_subset[cur_no_subsets+10]
        subset_extfaces = grid_ggd.grid_subset[cur_no_subsets+11]
        subset_extcells = grid_ggd.grid_subset[cur_no_subsets+12]
        subset_b25nodes.identifier.description = "All nodes from B2.5"
        subset_b25faces.identifier.description = "All faces from B2.5"
        subset_b25cells.identifier.description = "All cells from B2.5"
        subset_trinodes.identifier.index = -101
        subset_trinodes.identifier.name = "nodes_EIRENE"
        subset_trinodes.identifier.description = "Triangular mesh nodes from EIRENE"
        subset_trifaces.identifier.index = -102
        subset_trifaces.identifier.name = "faces_EIRENE"
        subset_trifaces.identifier.description = "Triangular mesh faces from EIRENE"
        subset_tricells.identifier.index = -105
        subset_tricells.identifier.name = "cells_EIRENE"
        subset_tricells.identifier.description = "Triangular mesh cells from EIRENE"
        subset_comnodes.identifier.index = -201
        subset_comnodes.identifier.name = "nodes_EIRENE_B2.5"
        subset_comnodes.identifier.description = "Triangular mesh nodes common between EIRENE and B2.5"
        subset_comfaces.identifier.index = -202
        subset_comfaces.identifier.name = "faces_EIRENE_B2.5"
        subset_comfaces.identifier.description = "Triangular mesh faces common between EIRENE and B2.5"
        subset_comcells.identifier.index = -205
        subset_comcells.identifier.name = "cells_EIRENE_B2.5"
        subset_comcells.identifier.description = "Triangular mesh cells overlapping between EIRENE and B2.5"
        subset_extnodes.identifier.index = -301
        subset_extnodes.identifier.name = "nodes_EIRENE_Extension_Only"
        subset_extnodes.identifier.description = "Triangular mesh nodes only in EIRENE extended mesh outside B2.5"
        subset_extfaces.identifier.index = -302
        subset_extfaces.identifier.name = "faces_EIRENE_Extension_Only"
        subset_extfaces.identifier.description = "Triangular mesh faces only in EIRENE extended mesh outside B2.5"
        subset_extcells.identifier.index = -305
        subset_extcells.identifier.name = "cells_EIRENE_Extension_Only"
        subset_extcells.identifier.description = "Triangular mesh cells only in EIRENE extended mesh outside B2.5"

        # Adding new node positions and cell corners data
        f33 = readdlm(fort[1])
        fnnodes = f33[1, 1]
        fnodeXnodeY = vec(f33[2:end, :]') * 1e-2 # cm to m
        fnodeX = fnodeXnodeY[1:fnnodes]
        fnodeY = fnodeXnodeY[fnnodes+1:end]
        fnode_inds = Array{Int64}(undef, fnnodes)
        for fnind ∈ 1:fnnodes
            i_existing = search_points(
                nodes,
                fnodeX[fnind],
                fnodeY[fnind];
                tol=fort_tol,
            )[1]
            if i_existing == 0
                resize!(nodes, length(nodes) + 1)
                this_node_ind = length(nodes)
                nodes[this_node_ind].geometry = [fnodeX[fnind], fnodeY[fnind]]

                add_subset_element!(subset_nodes, 1, 1, this_node_ind)
            else
                this_node_ind = i_existing[1]
            end
            fnode_inds[fnind] = this_node_ind
            add_subset_element!(subset_trinodes, 1, 1, this_node_ind)
        end
        f34 = readdlm(fort[2])
        fntria = f34[1, 1]
        fntriIndNodes = f34[2:end, :]
        fntri_inds = Array{Int64}(undef, fntria)
        for fntri ∈ 1:fntria
            resize!(cells, length(cells) + 1)
            this_cell_ind = length(cells)
            fntri_inds[fntri] = this_cell_ind
            resize!(cells[this_cell_ind].nodes, 3)
            cells[this_cell_ind].nodes = [
                fnode_inds[fntriIndNodes[fntri, 2]],
                fnode_inds[fntriIndNodes[fntri, 3]],
                fnode_inds[fntriIndNodes[fntri, 4]],
            ]
            resize!(cells[this_cell_ind].boundary, 3)
            for (bnd_ind, edge_pair) ∈ chosen_tri_edge_order
                tri_edge_nodes =
                    [cells[this_cell_ind].nodes[icorner] for icorner ∈ edge_pair]
                existing_edge_ind = search_edges(edges, tri_edge_nodes)
                if existing_edge_ind == 0
                    resize!(edges, length(edges) + 1)
                    this_edge_ind = length(edges)
                    edges[this_edge_ind].nodes = tri_edge_nodes
                    for (ii, edge_bnd) ∈ enumerate(edges[this_edge_ind].boundary)
                        edge_bnd.index = tri_edge_nodes[ii]
                    end
                    edges[this_edge_ind].measure =
                        distance_between_nodes(nodes, tri_edge_nodes)
                    add_subset_element!(subset_faces, 1, 2, this_edge_ind)
                else
                    this_edge_ind = existing_edge_ind
                end
                cells[this_cell_ind].boundary[bnd_ind].index = this_edge_ind
                add_subset_element!(subset_trifaces, 1, 2, this_edge_ind)
            end
            add_subset_element!(subset_cells, 1, 3, this_cell_ind)
            add_subset_element!(subset_tricells, 1, 3, this_cell_ind)
        end
        f35 = readdlm(fort[3])
        f35 = f35[2:end, :]
        for ii ∈ 1:fntria
            this_cell_ind = fntri_inds[f35[ii, 1]]
            this_cell = cells[this_cell_ind]
            for (bnd_no, bnd) ∈ enumerate(this_cell.boundary)
                neigh_col_ind = (bnd_no - 1) * 3 + 2
                if f35[ii, neigh_col_ind] != 0
                    bnd.neighbours = [fntri_inds[f35[ii, neigh_col_ind]]]
                else
                    bnd.neighbours = Int64[]
                end
            end
            if f35[ii, 11] != -1 && f35[ii, 12] != -1
                add_subset_element!(subset_comcells, 1, 3, this_cell_ind)
            else
                add_subset_element!(subset_extcells, 1, 3, this_cell_ind)
            end
        end

        subset_comnodes.element =
            subset_do(intersect, subset_b25nodes.element, subset_trinodes.element)
        subset_comfaces.element =
            subset_do(intersect, subset_b25faces.element, subset_trifaces.element)

        subset_extnodes.element =
            subset_do(setdiff, subset_trinodes.element, subset_b25nodes.element)
        subset_extfaces.element =
            subset_do(setdiff, subset_trifaces.element, subset_b25faces.element)
    end

    # Filling data in ggd now
    b2 = read_b2_output(b2output)
    # Add grid_ggd array equal to number of time steps
    resize!(ids.edge_profiles.ggd, b2["dim"]["time"])
    for it ∈ 1:b2["dim"]["time"]
        ggd = ids.edge_profiles.ggd[it]
        ggd.time = Float64.(b2["data"]["timesa"][it])
        grid_ggd_ind = gsdesc["identifier"]["index"]
        for (key, data) ∈ b2["data"]
            parent, prop = val_obj(ggd, key, grid_ggd_ind; gsi_ch=gsi_ch)
            if !isnothing(parent)
                setproperty!(parent, prop, data_xytoc(data[it, :, :]; nx=nx))
            end
        end
        # Done with filling data for this time step
    end # End of it

    # Adding magnetic field data
    if "bb" ∈ keys(gmtry["data"]) && load_bb
        bb = gmtry["data"]["bb"]
        if length(ids.equilibrium.time_slice) < gmtry["dim"]["time"]
            resize!(ids.equilibrium.time_slice, gmtry["dim"]["time"])
            ids.equilibrium.time = gmtry["data"]["timesa"]
        end
        for it ∈ 1:gmtry["dim"]["time"]
            # Note
            # Ideally, equilibrium keeps separate grids_ggd object for each time step
            # But since we have already created them in edge_profiles.grid_ggd, we
            # will not duplicate the information further.
            # If some other code requires it, it can done by
            # ids.equilibrium.grids_ggd = ids.edge_profiles.grid_ggd
            time_slice = ids.equilibrium.time_slice[it]
            resize!(time_slice.ggd, 1)
            resize!(time_slice.ggd[1].b_field_tor, 1)
            resize!(time_slice.ggd[1].b_field_r, 1)
            resize!(time_slice.ggd[1].b_field_z, 1)

            b_t = time_slice.ggd[1].b_field_tor[1]
            b_r = time_slice.ggd[1].b_field_r[1]
            b_z = time_slice.ggd[1].b_field_z[1]

            b_z.grid_index =
                b_r.grid_index = b_t.grid_index = gsdesc["identifier"]["index"]
            b_z.grid_subset_index = b_r.grid_subset_index = b_t.grid_subset_index = 5
            resize!(b_z.values, ncell)
            resize!(b_r.values, ncell)
            resize!(b_t.values, ncell)
            for iy ∈ 1:ny
                for ix ∈ 1:nx
                    ic::Int = (iy - 1) * nx + ix
                    b_z.values[ic] = bb[it, 1, iy, ix]
                    b_r.values[ic] = bb[it, 2, iy, ix]
                    b_t.values[ic] = bb[it, 3, iy, ix]
                end
            end
        end
    end
    return ids
end

end # module SOLPS2IMAS
