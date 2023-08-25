module SOLPS2IMAS

# using DelimitedFiles
using Revise
# using IMASDD
# const OMAS = IMASDD
using OMAS
using NCDatasets
using YAML

export try_omas
export populate_grid_ggd
export generate_test_data
export read_b2_output
export search_points
export solps2imas

function try_omas()
    dd = OMAS.dd()
    resize!(dd.equilibrium.time_slice, 1)
    dd.equilibrium.time_slice[1].profiles_1d.psi = [0.0, 1.0, 2.0, 3.9]
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

            if tag == "nx,ny,ns"  # This is present in b2fstate
                nx, ny, ns = array_line
                nx += 2  # Account for guard cells
                ny += 2  # Account for guard cells
                ret_dict["dim"] = Dict("nx" => nx, "ny" => ny, "ns" => ns)
                delete!(contents, "nx,ny,ns")
            elseif tag == "nx,ny"  # This is present in b2fgmtry
                nx, ny = array_line
                ns = 0
                nx += 2  # Account for guard cells
                ny += 2  # Account for guard cells
                ret_dict["dim"] = Dict("nx" => nx, "ny" => ny)
                delete!(contents, "nx,ny,ns")
            end
        end
    end
    # Cleanup arrays if applicable and return structured dictionary
    ret_dict["dim"]["time"] = 1  # This part of code will be for final state or geometry file
    ret_dict["data"] = Dict()    # Adding placeholder timestamp
    ret_dict["data"]["timesa"] = [0.0]
    for tag in keys(contents)
        # Note that size 1 dimention is added to the left always for time dimension
        if array_sizes[tag] == nx * ny
            ret_dict["data"][tag] = reshape(contents[tag], (1, ny, nx))
        elseif array_sizes[tag] == nx * ny * ns
            # If ns == 2, then r,z vector arrays can't be distinguished
            # from species-dependent quantities by their shapes. But
            # they get treated the same way, so it's okay.
            ret_dict["data"][tag] = reshape(contents[tag], (1, ns, ny, nx))
        elseif array_sizes[tag] == nx * ny * 2
            ret_dict["data"][tag] = reshape(contents[tag], (1, 2, ny, nx))
        elseif array_sizes[tag] == nx * ny * 2 * ns
            ret_dict["data"][tag] = reshape(contents[tag], (1, ns, 2, ny, nx))
        elseif array_sizes[tag] == nx * ny * 4
            # This case is only applicable to b2fgmtry, so ns will be 0 if this is
            # relevant.
            # Therefore, this case won't be inappropriately blocked by
            # ns * 2 when ns=2
            ret_dict["data"][tag] = reshape(contents[tag], (1, 4, ny, nx))
        elseif tag ∉ keys(ret_dict["dim"])
            ret_dict[tag] = contents[tag]
        end
    end
    return ret_dict
end

function search_points(dd, r, z)
    n = length(r)
    indices = zeros(Int, n)
    grid_number = 1
    space_number = 1
    subset_idx_node = 1
    # If an index remains at 0, it means the point in question was not found
    nodes = dd.edge_profiles.grid_ggd[grid_number].space[space_number].objects_per_dimension[subset_idx_node].object
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
solps_var_to_imas = YAML.load_file("$(@__DIR__)/solps_var_to_imas.yml")
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
        prop.grid_subset_index = path[gsi_ind]
        return path_to_obj(prop, prop_to_obj)
    end
end

"""
    solps2imas(b2gmtry, b2output, gsdesc)

Main function of the module. Takes in a geometry file and a
output file (either b2time or b2fstate) and a grid_ggd
description in the form of a Dict or filename to equivalent
YAML file. Returns data in OMAS.dd datastructure.
"""
function solps2imas(b2gmtry, b2output, gsdesc)
    gmtry = read_b2_output(b2gmtry)
    b2 =read_b2_output(b2output)
    nt = b2["dim"]["time"]
    times = b2["data"]["timesa"]
    nx = gmtry["dim"]["nx"]
    ny = gmtry["dim"]["ny"]
    crx = gmtry["data"]["crx"]
    cry = gmtry["data"]["cry"]
    ncell = nx * ny

    # Initialize an empty OMAS data structre
    dd = OMAS.dd()

    if typeof(gsdesc) == String
        gsdesc = YAML.load_file(gsdesc)
    end

    # Add ggd and grid_ggd array equal to number of time steps
    resize!(dd.edge_profiles.ggd, nt)
    resize!(dd.edge_profiles.grid_ggd, nt)
    for it in 1:nt
        # Setup the grid first for this time step
        grid_ggd = dd.edge_profiles.grid_ggd[it]
        grid_ggd.time = Float64.(times[it])
        dict2prop(grid_ggd, gsdesc)
        for sn in keys(gsdesc["space"])
            space = grid_ggd.space[sn]
            # Assuming following to be standard for now. We can add this info through YAML as well
            resize!(space.objects_per_dimension, 4)
            o0 = space.objects_per_dimension[1]  # 0D objects
            o1 = space.objects_per_dimension[2]  # 1D objects
            o2 = space.objects_per_dimension[3]  # 2D objects
            o3 = space.objects_per_dimension[4]  # 3D objects

            # Resizing objects to hold cell geometry data
            # Should be fewer than this many points, but this way we won't under-fill
            resize!(o0.object, ncell * 4)  # Points
            resize!(o1.object, ncell * 4)  # Edges
            resize!(o2.object, ncell)  # Faces
            resize!(o3.object, ncell)  # Volumes

            # Initialize geometry for 0D objects
            for i = 1:(ncell * 4)
                o0.object[i].geometry = [0.0, 0.0]
            end
            # Initialize nodes for cells
            for i = 1:(ncell)
                o2.object[i].nodes = [0, 0, 0, 0]
            end

            j = 1
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
                        i_existing = search_points(dd, crx[1, icorner, iy, ix], cry[1, icorner, iy, ix])[1]
                        if i_existing == 0
                            o0.object[j].geometry = [crx[1, icorner, iy, ix], cry[1, icorner, iy, ix]]
                            o2.object[ic].nodes[icorner] = j
                            j += 1
                        else
                            o2.object[ic].nodes[icorner] = i_existing[1]
                        end
                    end
                end
            end
        end  # End of setting up space

        # Filling data in ggd now
        ggd = dd.edge_profiles.ggd[it]
        ggd.time = Float64.(times[it])
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
    # Store magnetic field data in equilibrium
    equilibrium_slice = 1  # This will be the case unless someone is doing something funny
    #if length(dd.equilibrium.time_slice) < equilibrium_slice:
    resize!(dd.equilibrium.time_slice, equilibrium_slice)
    #end
    #if length(dd.equilibrium.time_slice[equilibrium_slice].ggd) < 1
    resize!(dd.equilibrium.time_slice[equilibrium_slice].ggd, 1)
    #end
    eq_ggd = dd.equilibrium.time_slice[equilibrium_slice].ggd[1]
    resize!(eq_ggd.b_field_tor, 5)
    resize!(eq_ggd.b_field_tor[5].values, ncell)
    eq_bb = eq_ggd.b_field_tor[5].values
    data = gmtry["data"]["bb"][1, 3, :, :]  # bb stores Bpoloidal, Bradial, Btoroidal, Btotal;
    # Bpol and Brad are w.r.t. the mesh so Brad is always 0. It's not B_R.
    for iy = 1:ny
        for ix = 1:nx
            ic::Int = (iy - 1) * nx + ix
            eq_bb[ic] = data[iy, ix]
        end
    end
    return dd
end

end # module SOLPS2IMAS
