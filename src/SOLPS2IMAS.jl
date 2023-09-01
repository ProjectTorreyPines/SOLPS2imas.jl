module SOLPS2IMAS

using Revise
using OMAS: dd as OMAS_dd
using NCDatasets: Dataset, dimnames
using YAML: load_file as YAML_load_file

export try_omas
export generate_test_data
export read_b2_output
export search_points
export solps2imas


function try_omas()
    ids = OMAS_dd()
    resize!(ids.equilibrium.time_slice, 1)
    ids.equilibrium.time_slice[1].profiles_1d.psi = [0.0, 1.0, 2.0, 3.9]
    return nothing
end


function generate_test_data(nx=94, ny=38, sep=15, lcut=25, rcut=69)
    """This doesn't work well"""
    center_r = 1.7
    center_z = 0.0
    aminor = 0.6
    inner_minor = 0.3
    outer_minor = 0.8
    elongation = 1.8
    θₓ = 8 * π/6
    θᵢ = θₓ - π/4
    θₒ = θₓ + π/4
    inner_leg_length = 0.5
    outer_leg_length = 0.6
    dx_il = inner_leg_length / lcut
    dx_ol = outer_leg_length / (nx - rcut)
    dr_in = (aminor - inner_minor) / sep
    dr_out = (outer_minor = aminor) / (ny - sep)
    cell_centers_r = Array{Float64}(undef, ny, nx)
    cell_centers_z = Array{Float64}(undef, ny, nx)
    crx = Array{Float64}(undef, 4, ny, nx)
    cry = Array{Float64}(undef, 4, ny, nx)
    xpoint_r = center_r + aminor * cos(θₓ)
    xpoint_z = center_z + aminor * sin(θₓ) * elongation
    for ix::Int = 1:nx
        for iy::Int = 1:ny
            if iy < sep
                dy = dr_in
            else
                dy = dr_out
            end

            if (ix >= lcut) && (ix <= rcut)
                a = aminor + (iy - sep + 0.5) * dy
                θ = θₓ - (ix - lcut) / (rcut - lcut) * 2 * π
                # println(θ)
                cell_centers_r[iy, ix] = cos(θ) * a + center_r
                cell_centers_z[iy, ix] = sin(θ) * a * elongation + center_z
            elseif ix < lcut
                cell_centers_r[iy, ix] = xpoint_r + cos(θᵢ) * dx_il * (lcut - ix) + sin(θᵢ) * dy * (iy - sep)
                cell_centers_z[iy, ix] = xpoint_z + sin(θᵢ) * dx_il * (lcut - ix) + cos(θᵢ) * dy * (iy - sep)
            else
                cell_centers_r[iy, ix] = xpoint_r + cos(θₒ) * dx_ol * (ix - rcut) + sin(θₒ) * dy * (iy - sep)
                cell_centers_z[iy, ix] = xpoint_z + sin(θₒ) * dx_ol * (ix - rcut) + cos(θₒ) * dy * (iy - sep)
            end
        end
    end
    # Sanitize
    rmin = 0.0
    rmax = 3.0
    zmin = -3.0
    zmax = 3.0
    for i::Int = 1:nx*ny
        # if mod(i, nx) == 0
        #     print(cell_centers_r[i], " ", cell_centers_r[i] > rmax, " ", rmax)
        # end
        cell_centers_r[i] = cell_centers_r[i] > rmax ? rmax : cell_centers_r[i]
        cell_centers_z[i] = cell_centers_z[i] > zmax ? zmax : cell_centers_z[i]
        cell_centers_r[i] = cell_centers_r[i] < rmin ? rmin : cell_centers_r[i]
        cell_centers_z[i] = cell_centers_z[i] < zmin ? zmin : cell_centers_z[i]
        # if (mod(i, nx) == 0)
        #     println(" ", cell_centers_r[i])
        # end
    end
    return (cell_centers_r, cell_centers_z)
end


function read_b2time_output(filename)
    dim_order = (
        "time",
        "ns",
        "nstrat",
        "nc",
        "ndir",
        "ny", "nybl", "nybr", "nya", "nyi",
        "nx", "nxbl", "nxbr", "nxa", "nxi",
    )
    ret_dict = Dict("dim" => Dict(), "data" => Dict())
    ds = Dataset(filename)
    for key in keys(ds.dim)
        ret_dict["dim"][key] = ds.dim[key]
    end
    for key in keys(ds)
        if key != "ntstep"
            d = dimnames(ds[key])
            permute = [y for y in [findfirst(x->x==dimord, d) for dimord in dim_order] if y !== nothing]
            try
                ret_dict["data"][key] = permutedims(Array(ds[key]), permute)
            catch e
                println("Error in reading ", key)
                showerror(stdout, e)
                println("Continuing by ignoring this field")
            end
        end
    end
    return ret_dict
end


function read_b2_output(filename)
    # x = readdlm(filename, ' ', Float64, '\n', header=true)  # Doesn't handle the text lines at the start of each array

    if cmp(splitext(filename)[2], ".nc") == 0
        return read_b2time_output(filename)
    end

    contents = Dict()
    array_sizes = Dict()
    ret_dict = Dict()
    lines = open(filename) do f
        readlines(f)
    end
    nx = 0
    ny = 0
    ns = 0
    tag = ""
    arraysize = 0
    arraytype = nothing
    j = 1
    for l in lines
        if startswith(l, "*cf:")
            j = 1  # Reset intra-array element counter
            _, arraytype, arraysize, tag = split(l)
            arraysize = parse(Int, arraysize)
            if arraytype == "char"
                contents[tag] = ""
            elseif arraytype == "int"
                contents[tag] = Array{Int}(undef, arraysize)
            else
                contents[tag] = Array{Float64}(undef, arraysize)
            end
            array_sizes[tag] = arraysize
        elseif tag != ""
            if arraytype == "int"
                array_line = [parse(Int, ss) for ss in split(l)]
                array_inc = size(array_line)[1]
            elseif arraytype == "real"
                array_line = [parse(Float64, ss) for ss in split(l)]
                array_inc = size(array_line)[1]
            else
                array_line = l
                array_inc = 1
            end
            if arraytype == "char"
                contents[tag] = array_line
            else
                contents[tag][j:j+array_inc-1] = array_line
            end
            j += array_inc
        end
    end
    if "nx,ny" ∈ keys(contents)
        return extract_geometry(contents)
    elseif "nx,ny,ns" ∈ keys(contents)
        return extract_state_quantities(contents)
    else
        throw(DomainError(keys(contents),
              "nx,ny (b2fgmtry) or nx,ny,ns (b2fstate) must be present in b2 output file."))
    end
end


function extract_geometry(gmtry)
    ret_dict = Dict("dim" => Dict(), "data" => Dict())
    ret_dict["dim"]["nx_no_guard"], ret_dict["dim"]["ny_no_guard"] = gmtry["nx,ny"]
    # includes guard cells
    nx = ret_dict["dim"]["nx"] = ret_dict["dim"]["nx_no_guard"] + 2
    ny = ret_dict["dim"]["ny"] = ret_dict["dim"]["ny_no_guard"] + 2
    if "nncut" ∈ keys(gmtry)
        ret_dict["dim"]["nncut"] = gmtry["nncut"]
    end
    # Adding placeholder timestamp
    ret_dict["dim"]["time"] = 1
    ret_dict["data"]["timesa"] = [0.0]
    for k in keys(gmtry)
        # The 4 fields of bb are poloidal, radial, toroidal, and total magnetic field
        # according to page 212 of D. Coster, "SOLPS-ITER [manual]" (2019)
        # The 4 fields in crx and cry are the corners of each grid cell.
        if k ∈ ["crx", "cry", "bb"]
            ret_dict["data"][k] = permutedims(reshape(gmtry[k], (nx, ny, 4, 1)), (4, 3, 2, 1))
        elseif k ∈ ["leftcut", "bottomcut", "rightcut", "topcut"]
            ret_dict["data"][k] = Array([gmtry[k][1]])
        elseif k ∈ ["leftcut2", "bottomcut2", "rightcut2", "topcut2"]
            ret_dict["data"][k] = Array([gmtry[k[1:end-1]][1], gmtry[k][1]])
        elseif length(gmtry[k]) == nx * ny
            ret_dict["data"][k] = permutedims(reshape(gmtry[k], (nx, ny, 1)), (3, 2, 1))
        elseif k ∉ keys(ret_dict["dim"])
            ret_dict["data"][k] = gmtry[k]
        end
    end
    return ret_dict
end


function extract_state_quantities(state)
    ret_dict = Dict("dim" => Dict(), "data" => Dict())
    ret_dict["dim"]["nx_no_guard"], ret_dict["dim"]["ny_no_guard"], ret_dict["dim"]["ns"] = state["nx,ny,ns"]
    # includes guard cells
    nx = ret_dict["dim"]["nx"] = ret_dict["dim"]["nx_no_guard"] + 2
    ny = ret_dict["dim"]["ny"] = ret_dict["dim"]["ny_no_guard"] + 2
    ns = ret_dict["dim"]["ns"]
    ndir = ret_dict["dim"]["ndir"] = 2
    # Adding placeholder timestamp
    ret_dict["dim"]["time"] = 1
    ret_dict["data"]["timesa"] = [0.0]
    for k in keys(state)
        l = length(state[k])
        if l == nx * ny
            ret_dict["data"][k] = permutedims(reshape(state[k], (nx, ny, 1)), (3, 2, 1))
        elseif l == nx * ny * ns
            ret_dict["data"][k] = permutedims(reshape(state[k], (nx, ny, ns, 1)), (4, 3, 2, 1))
        elseif l == nx * ny * ndir
            ret_dict["data"][k] = permutedims(reshape(state[k], (nx, ny, ndir, 1)), (4, 3, 2, 1))
        elseif l == nx * ny * ndir * ns
            ret_dict["data"][k] = permutedims(reshape(state[k], (nx, ny, ndir, ns, 1)), (5, 4, 3, 2, 1))
        elseif l == ns
            ret_dict["data"][k] = permutedims(reshape(state[k], (ns, 1)), (2, 1))
        elseif k ∉ keys(ret_dict["dim"])
            ret_dict["data"][k] = state[k]
        end
    end
    return ret_dict
end


"""
   select_core(q, topcut, bottomcut, leftcut, rightcut)

Selects core part for any quantity q
Inspired from
https://odin.gat.com/eldond/sparc_detach_ctrl/blob/master/synth_diag/scripts/profile_extensions.py#L20
"""
function select_core(q; topcut, bottomcut, leftcut, rightcut)
    slices = [range(1, s) for s in size(q)[1:end-2]]
    append!(slices, [range(bottomcut + 2, topcut + 1), range((leftcut + 1), (rightcut + 1))])
    return view(q, slices...)
end


"""
   get_cell_numbers(nx, ny; topcut, bottomcut, leftcut, rightcut)

Returns number of cells in core, SOL, inner divertor region and outer divertor region
"""
function get_cell_numbers(nx, ny; topcut, bottomcut, leftcut, rightcut)
    ncell_core = (topcut - bottomcut) * (rightcut - leftcut + 1)
    ncell_sol = (ny - topcut - 1) * nx
    ncell_idr = (topcut - bottomcut) * leftcut
    ncell_odr = (topcut - bottomcut) * (nx - rightcut - 1)
    return ncell_core, ncell_sol, ncell_idr, ncell_odr
end


"""
   in_core(ix, iy; topcut, bottomcut, leftcut, rightcut)

Returns true if cell indexed ix, iy lie inside the core
"""
function in_core(ix, iy; topcut, bottomcut, leftcut, rightcut)
    return bottomcut + 1 < iy < topcut + 2 && leftcut < ix < rightcut + 2
end


"""
   in_sol(ix, iy; topcut, bottomcut, leftcut, rightcut)

Returns true if cell indexed ix, iy lie inside the SOL
"""
function in_sol(ix, iy; topcut, bottomcut, leftcut, rightcut)
    return topcut + 1 < iy
end


"""
   in_idr(ix, iy; topcut, bottomcut, leftcut, rightcut)

Returns true if cell indexed ix, iy lie inside the inner divertor region
"""
function in_idr(ix, iy; topcut, bottomcut, leftcut, rightcut)
    return bottomcut + 1 < iy < topcut + 2 && ix < leftcut + 1
end


"""
   in_odr(ix, iy; topcut, bottomcut, leftcut, rightcut)

Returns true if cell indexed ix, iy lie inside the outer divertor region
"""
function in_odr(ix, iy; topcut, bottomcut, leftcut, rightcut)
    return bottomcut + 1 < iy < topcut + 2 && rightcut + 1 < ix
end


"""
   add_subset_element!(subset, dd_ind, sn, dim, index, ix, iy, in_subset=(x...)->true; kwargs...)

Adds the geometric element in subset object (assumed to be resized already) at element dd index dd_ind,
with space number sn, dimension dim, index index, for a cell in SOLPS mesh indices ix, iy. To determine,
if the element should be added or not, a function in_subset can be provided that gets the arguments
(ix, iy; kwargs...). These functions will be in_core, in_sol etc as difined above.
"""
function add_subset_element!(subset, dd_ind, sn, dim, index, ix, iy, in_subset=(x...)->true; kwargs...)
    if in_subset(ix, iy; kwargs...)
        dd_ind += 1
        resize!(subset.element[dd_ind].object, 1)
        subset.element[dd_ind].object[1].space = sn
        subset.element[dd_ind].object[1].dimension = dim
        subset.element[dd_ind].object[1].index = index
    end
end


function search_points(ids, r, z)
    n = length(r)
    indices = zeros(Int, n)
    grid_number = 1
    space_number = 1
    subset_idx_node = 1
    # If an index remains at 0, it means the point in question was not found
    nodes = ids.edge_profiles.grid_ggd[grid_number].space[space_number].objects_per_dimension[subset_idx_node].object
    for j = 1:n
        for i in eachindex(nodes)
            rn = nodes[i].geometry[1]
            zn = nodes[i].geometry[2]
            if (rn == r[j]) && (zn == z[j])
                indices[j] = i
                break
            end
        end
    end
    return indices
end


"""
    dict2prop(obj, dict)

Copies grid_ggd and space description in dict format to the data structure recursively.
"""
function dict2prop(obj, dict)
    for (key, prop) in dict
        if isa(key, Int)
            resize!(obj, key)
            dict2prop(obj[key], prop)
        elseif isa(key, String)
            if isa(prop, Dict)
                dict2prop(getfield(obj, Symbol(key)), prop)
            else
                setproperty!(obj, Symbol(key), prop)
            end
        end
    end
end


"""
    path_to_obj(obj, path)

Using path which is a collection of strings and integers,
return a field of structure obj, using strings as names
of fields and integers as index of an array of structure.
Example:
path_to_obj(obj, ["abc", 3, "cde", "fgh", 5, "ijk"])
returns
obj.abc[3].cde.fgh[5].ijk
Note: If integer is -1, the array of field is resized to
increase by 1 in length and last element is returned.
"""
function path_to_obj(obj, path)
    for ele in path
        if typeof(ele) == String
            obj = getfield(obj, Symbol(ele))
        elseif typeof(ele) == Int
            if ele == -1
                resize!(obj, length(obj) + 1)
                ele = length(obj)
            elseif length(obj) < ele
                resize!(obj, ele)
            end
            obj = obj[ele]
        end
    end
    return obj
end


isint(x) = typeof(x) == Int


"""
    val_obj(var, ggd, grid_ggd_index)

Given SOLPS variable name (var), return the field to write values on from
ggd object. 
"""
solps_var_to_imas = YAML_load_file("$(@__DIR__)/solps_var_to_imas.yml")
function val_obj(ggd, var, grid_ggd_index)
    if var ∉ keys(solps_var_to_imas)
        return nothing
    else
        path = solps_var_to_imas[var]
        gsi_ind = findlast(isint, path)
        path_to_prop = path[1:gsi_ind]
        prop_to_obj = path[gsi_ind + 1:end]
        prop = path_to_obj(ggd, path_to_prop)
        prop.grid_index = grid_ggd_index
        prop.grid_subset_index = gsi_ind
        return path_to_obj(prop, prop_to_obj)
    end
end


"""
    find_subset_index()

Finds the julia index of the subset with ggd_index matching the request.
Example: GGD defines subset index 5 as being all 2D cells. But what is the Julia
index of that subset within the IMAS DD representation? Yes, we're trying to find
the index of the index, but they're different meanings of index.
We'll call them dd_index (the place in the DD) and ggd_index (the type of subset).
"""
function find_subset_index(gsdesc, ggd_index)
    subsets = length(gsdesc["grid_subset"])
    for dd_index = 1:subsets
        if gsdesc["grid_subset"][dd_index]["identifier"]["index"] == ggd_index
            return dd_index
        end
    end
    return 0  # Indicates failure
end


"""
    solps2imas(b2gmtry, b2output, gsdesc)

Main function of the module. Takes in a geometry file and a
output file (either b2time or b2fstate) and a grid_ggd
description in the form of a Dict or filename to equivalent
YAML file. Returns data in OMAS.dd datastructure.
"""
function solps2imas(b2gmtry, b2output, gsdesc; load_bb=false)
    # Initialize an empty OMAS data structre
    ids = OMAS_dd()

    # Setup the grid first
    gmtry = read_b2_output(b2gmtry)

    nx = gmtry["dim"]["nx"]
    ny = gmtry["dim"]["ny"]
    ncell = nx * ny
    crx = gmtry["data"]["crx"]
    cry = gmtry["data"]["cry"]
    cut_keys = ["leftcut", "rightcut", "bottomcut", "topcut"]
    cuts_found = cut_keys ⊆ keys(gmtry["data"])
    if cuts_found
        cuts = Dict([(Symbol(key), gmtry["data"][key][1]) for key in cut_keys])
        ncell_core, ncell_sol, ncell_idr, ncell_odr = get_cell_numbers(nx, ny; cuts...)
    end


    if typeof(gsdesc) == String
        gsdesc = YAML_load_file(gsdesc)
    end

    # Add grid_ggd array equal to number of time steps
    resize!(ids.edge_profiles.grid_ggd, gmtry["dim"]["time"])
    for it in 1:gmtry["dim"]["time"]
        
        grid_ggd = ids.edge_profiles.grid_ggd[it]
        grid_ggd.time = Float64.(gmtry["data"]["timesa"][it])
        dict2prop(grid_ggd, gsdesc)
        for sn in keys(gsdesc["space"])
            space = grid_ggd.space[sn]
            # Assuming following to be standard for now. We can add this info through YAML as well
            resize!(space.objects_per_dimension, 4)
            o0 = space.objects_per_dimension[1]  # 0D objects
            o1 = space.objects_per_dimension[2]  # 1D objects
            o2 = space.objects_per_dimension[3]  # 2D objects
            o3 = space.objects_per_dimension[4]  # 3D objects

            subsets = length(gsdesc["grid_subset"])
            resize!(grid_ggd.grid_subset, subsets)
            for isub = 1:subsets
                grid_ggd.grid_subset[isub].identifier.name = gsdesc["grid_subset"][isub]["identifier"]["name"]
                grid_ggd.grid_subset[isub].identifier.index = gsdesc["grid_subset"][isub]["identifier"]["index"]
                grid_ggd.grid_subset[isub].identifier.description = gsdesc["grid_subset"][isub]["identifier"]["description"]
                grid_ggd.grid_subset[isub].dimension = gsdesc["grid_subset"][isub]["dimension"]
            end

            subset_nodes = grid_ggd.grid_subset[find_subset_index(gsdesc, 1)]  # nodes have index 1
            subset_faces = grid_ggd.grid_subset[find_subset_index(gsdesc, 2)]  # faces (edges with elongation in third dimension) have index 2
            subset_cells = grid_ggd.grid_subset[find_subset_index(gsdesc, 5)]  # cells (cell in 2D grid, volume with elongation) have index 5
            if cuts_found
                subset_core = grid_ggd.grid_subset[find_subset_index(gsdesc, 22)]  # core cells have index 22
                subset_sol = grid_ggd.grid_subset[find_subset_index(gsdesc, 23)]   # sol cells have index 23
                subset_odr = grid_ggd.grid_subset[find_subset_index(gsdesc, 24)]   # odr cells have index 24
                subset_idr = grid_ggd.grid_subset[find_subset_index(gsdesc, 25)]   # idr cells have index 25
            end

            # Resizing objects to hold cell geometry data
            # Should be fewer than this many points, but this way we won't under-fill
            resize!(o0.object, ncell * 4)  # Points
            resize!(subset_nodes.element, ncell * 4)
            resize!(o1.object, ncell * 4)  # Faces / edges
            resize!(subset_faces.element, ncell * 4)
            resize!(o2.object, ncell)  # Cells (2D)
            resize!(subset_cells.element, ncell)
            resize!(o3.object, ncell)  # Volumes
            if cuts_found
                resize!(subset_core.element, ncell_core)
                resize!(subset_sol.element, ncell_sol)
                resize!(subset_idr.element, ncell_idr)
                resize!(subset_odr.element, ncell_odr)
            end

            # Initialize geometry for 0D objects
            for i = 1:(ncell * 4)
                o0.object[i].geometry = [0.0, 0.0]
            end
            # Initialize nodes for cells
            for i = 1:(ncell)
                o2.object[i].nodes = [0, 0, 0, 0]
            end

            j = 1
            cell_dd_ind = 0
            core_dd_ind = 0
            sol_dd_ind = 0
            idr_dd_ind = 0
            odr_dd_ind = 0
            for iy = 1:ny
                for ix = 1:nx
                    ic::Int = (iy - 1) * nx + ix
                    # Adding node positions data to grid_ggd[grid_number].space[space_number].objects_per_dimension[0].object[:].geometry
                    # Adding cell corners data to grid_ggd[grid_number].space[space_number].objects_per_dimension[2].object[:].nodes[1:4]
                    for icorner = 1:4
                        # Have to search to see if the node is already added and then record its index
                        # If not already listed, then list it under new index and record that
                        # Note that time index has been fixed to 1 here. Only handling fixed grid geometry
                        # through the run cases.
                        i_existing = search_points(ids, crx[1, icorner, iy, ix], cry[1, icorner, iy, ix])[1]
                        if i_existing == 0
                            o0.object[j].geometry = [crx[1, icorner, iy, ix], cry[1, icorner, iy, ix]]
                            o2.object[ic].nodes[icorner] = j
                            resize!(subset_nodes.element[j].object, 1)
                            subset_nodes.element[j].object[1].space = sn
                            subset_nodes.element[j].object[1].dimension = 0
                            subset_nodes.element[j].object[1].index = j
                            j += 1
                        else
                            o2.object[ic].nodes[icorner] = i_existing[1]
                        end
                    end
                    add_subset_element!(subset_cells, cell_dd_ind, sn, 2, ic, ix, iy)
                    if cuts_found
                        # add_subset_element!(subset, dd_ind, sn, dim, index, ix, iy, in_subset=(x...) -> true; kwargs...)
                        add_subset_element!(subset_core, core_dd_ind, sn, 2, ic, ix, iy, in_core; cuts...)
                        add_subset_element!(subset_sol, sol_dd_ind, sn, 2, ic, ix, iy, in_sol; cuts...)
                        add_subset_element!(subset_idr, idr_dd_ind, sn, 2, ic, ix, iy, in_idr; cuts...)
                        add_subset_element!(subset_odr, odr_dd_ind, sn, 2, ic, ix, iy, in_odr; cuts...)
                    end
                end
            end
        end  # End of setting up space
    end

    # Filling data in ggd now
    b2 = read_b2_output(b2output)
    # Add grid_ggd array equal to number of time steps
    resize!(ids.edge_profiles.ggd, b2["dim"]["time"])
    for it in 1:b2["dim"]["time"]
        ggd = ids.edge_profiles.ggd[it]
        ggd.time = Float64.(b2["data"]["timesa"][it])
        for (key, data) in b2["data"]
            obj = val_obj(ggd, key, gsdesc["identifier"]["index"])
            if !isnothing(obj)
                resize!(obj, ncell)
                for iy = 1:ny
                    for ix = 1:nx
                        ic::Int = (iy - 1) * nx + ix
                        obj[ic] = data[it, iy, ix]
                    end
                end
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
        for it in 1:gmtry["dim"]["time"]
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

            b_z.grid_index = b_r.grid_index = b_t.grid_index = gsdesc["identifier"]["index"]
            b_z.grid_subset_index = b_r.grid_subset_index = b_t.grid_subset_index = 5
            resize!(b_z.values, ncell)
            resize!(b_r.values, ncell)
            resize!(b_t.values, ncell)
            for iy = 1:ny
                for ix = 1:nx
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
