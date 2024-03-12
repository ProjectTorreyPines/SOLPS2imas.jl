module SOLPS2IMAS

using Revise
using IMASDD: IMASDD
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

"""
    val_obj(var, ggd, grid_ggd_index)

Given SOLPS variable name (var), returns pair of parent object and property name
to write value on. If var is not found in solps_var_to_imas, returns nothing, nothing.
"""
solps_var_to_imas = YAML_load_file("$(@__DIR__)/solps_var_to_imas.yml")
function val_obj(ggd, var, grid_ggd_index)
    if var ∉ keys(solps_var_to_imas)
        return nothing, nothing
    else
        path, gsi = solps_var_to_imas[var]
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

"""
    solps2imas(b2gmtry, b2output, gsdesc; load_bb=false)

Main function of the module. Takes in a geometry file and a
output file (either b2time or b2fstate) and a grid_ggd
description in the form of a Dict or filename to equivalent
YAML file. Returns data in IMASDD.dd datastructure.
"""
function solps2imas(
    b2gmtry,
    b2output;
    gsdesc=default_gsdesc,
    b2mn=nothing,
    load_bb=false,
)
    # Initialize an empty IMAS data structre
    ids = IMASDD.dd()

    # Setup the grid first
    gmtry = read_b2_output(b2gmtry)

    jxi = jxa = nothing
    if !isnothing(b2mn)
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
            nodes = resize!(o1.object, ncell * 4)  # Nodes (1D)
            edges = resize!(o2.object, ncell * 4)  # Edges (2D)
            cells = resize!(o3.object, ncell)  # Cells (3D)

            # Initialize geometry for 1D objects(nodes), nodes for 2D objects(edges)
            for i ∈ 1:(ncell*4)
                nodes[i].geometry = [0.0, 0.0]
                edges[i].nodes = [0, 0]
                resize!(edges[i].boundary, 2)
                for bnd ∈ edges[i].boundary
                    bnd.neighbours = Int64[]
                end
            end
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
                            edges[edge_ind].nodes = edge_nodes
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

    # Filling data in ggd now
    b2 = read_b2_output(b2output)
    # Add grid_ggd array equal to number of time steps
    resize!(ids.edge_profiles.ggd, b2["dim"]["time"])
    for it ∈ 1:b2["dim"]["time"]
        ggd = ids.edge_profiles.ggd[it]
        ggd.time = Float64.(b2["data"]["timesa"][it])
        for (key, data) ∈ b2["data"]
            parent, prop = val_obj(ggd, key, gsdesc["identifier"]["index"])
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
