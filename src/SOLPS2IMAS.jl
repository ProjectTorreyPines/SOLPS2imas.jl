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
        "ndir",
        "ny", "nybl", "nybr", "nya", "nyi",
        "nx", "nxbl", "nxbr", "nxa", "nxi",
    )
    contents = Dict()
    ds = Dataset(filename)
    for key in keys(ds.dim)
        contents[key] = ds.dim[key]
    end
    for key in keys(ds)
        contents[key] = Array(ds[key])
        d = dimnames(ds[key])
        permute = [y for y in [findfirst(x->x==dimord, d) for dimord in dim_order] if y !== nothing]
        #println(key)
        #println(d, size(contents[key]))
        #println(permute)
        #print(key, " ", size(contents[key]), " ")
        contents[key] = permutedims(contents[key], permute)
        #println(size(contents[key]))
    end
    return contents
end

function read_b2_output(filename)
    # x = readdlm(filename, ' ', Float64, '\n', header=true)  # Doesn't handle the text lines at the start of each array

    if cmp(splitext(filename)[2], ".nc") == 0
        return read_b2time_output(filename)
    end

    contents = Dict()
    array_sizes = Dict()
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
            elseif tag == "nx,ny"  # This is present in b2fgmtry
                nx, ny = array_line
                ns = 0
                nx += 2  # Account for guard cells
                ny += 2  # Account for guard cells
            end
        end
    end
    # Cleanup arrays if applicable
    for tag in keys(contents)
        if array_sizes[tag] == nx * ny
            # println(tag, " ", array_sizes[tag], " ", nx, " ", ny, " ", size(contents[tag]))
            # print(contents[tag])
            contents[tag] = reshape(contents[tag], (ny, nx))
        elseif array_sizes[tag] == nx * ny * ns
            # If ns == 2, then r,z vector arrays can't be distinguished
            # from species-dependent quantities by their shapes. But
            # they get treated the same way, so it's okay.
            contents[tag] = reshape(contents[tag], (ns, ny, nx))
        elseif array_sizes[tag] == nx * ny * 2
            contents[tag] = reshape(contents[tag], (2, ny, nx))
        elseif array_sizes[tag] == nx * ny * 2 * ns
            contents[tag] = reshape(contents[tag], (ns, 2, ny, nx))
        elseif array_sizes[tag] == nx * ny * 4
            # This case is only applicable to b2fgmtry, so ns will be 0 if this is
            # relevant.
            # Therefore, this case won't be inappropriately blocked by
            # ns * 2 when ns=2
            contents[tag] = reshape(contents[tag], (4, ny, nx))
        end
    end
    return contents
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
        for i = 1:length(nodes)
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
    populate_grid_ggd(nx::Int64, ny::Int64, crx, cry, group, quantity, values, times)

Create a OMAS.dd datastructure from given SOLPS data for one quantity of interest.
Arguments:
gsdesc: Dict or String with filename to YAML file containing grid_ggd and space description
crx: Array(4, nx, ny) of r-coordinates of cell corners (cornerIndex x cellXIndex x cellYIndex)
cry: Array(4, nx, ny) of z-coordinates of cell corners
group: Field group (e_field, electrons, ion, neutral etc)
quantity: Data quantity specifier (density, velocity, pressure etc)
values: Array(len(times), nx, ny) of values of the quantity in a cell at a time
times: Array of times if data contains third dimension for time evolution
dd: (For future, Optional) OMAS.dd datastructure to add data to. Creates empty new structure if not provided.

Returns:
dd: OMAS.dd Data structure
"""
function populate_grid_ggd(gsdesc, crx, cry, group, quantity, values, times)
    _, ny, nx = size(crx)
    ncell = nx * ny
    dd = OMAS.dd()
    ndim = length(size(values))
    if times !== nothing
        nt = length(times)
        time_data_present = true
    else
        nt = 1
        time_data_present = false
    end

    if typeof(gsdesc) == String
        gsdesc = YAML.load_file(gsdesc)
    end
    for grid_number in keys(gsdesc)
        # Add new grid_ggd at the specified grid_number
        if length(dd.edge_profiles.grid_ggd) < grid_number
            resize!(dd.edge_profiles.grid_ggd, grid_number)
        end

        grid_ggdn = dd.edge_profiles.grid_ggd[grid_number]
        for pn in propertynames(grid_ggdn.identifier)
            if pn in keys(gsdesc[grid_number])
                setproperty!(grid_ggdn.identifier, pn, gsdesc[grid_number][pn])
            end
        end

        # id = grid_ggdn.identifier
        # id.name = gsdesc[grid_number]["name"]
        # id.index = gsdesc[grid_number]["index"]
        # id.description = gsdesc[grid_number]["description"]

        for space_number in keys(gsdesc[grid_number]["space"])
            # Add new space at the specified space_number
            if length(grid_ggdn.space) < space_number
                resize!(grid_ggdn.space, space_number)
            end

            space = grid_ggdn.space[space_number]

            for pn in propertynames(space.identifier)
                if pn in keys(gsdesc[grid_number]["space"][space_number])
                    setproperty!(space.identifier, pn, gsdesc[grid_number]["space"][space_number][pn])
                end
            end

            for pn in propertynames(space.geometry_type)
                if pn in keys(gsdesc[grid_number]["space"][space_number]["geometry_type"])
                    setproperty!(space.identifier, pn, gsdesc[grid_number]["space"][space_number]["geometry_type"][pn])
                end
            end

            space.coordinates_type = gsdesc[grid_number]["space"][space_number]["coordinates_type"]

            # space.identifier.name = "sp4ce"
            # space.identifier.index = 1
            # space.identifier.description = "The final frontier"
            # space.geometry_type.name = "standard"  # I doubt this is needed
            # space.geometry_type.index = 0  # 0 for standard, 1 for fourier. This is the important field
            # space.geometry_type.description = "trying to hold a b2/solps mesh here"  # I doubt this is needed
            # space.coordinates_type = [4, 3]  # r, z

            # Assuming following to be standard for now. We can add this info through YAML as well
            resize!(space.objects_per_dimension, 4)
            o0 = space.objects_per_dimension[1]  # 0D objects
            o1 = space.objects_per_dimension[2]  # 1D objects
            o2 = space.objects_per_dimension[3]  # 2D objects
            o3 = space.objects_per_dimension[4]  # 3D objects
            subset_idx_node = 1
            subset_idx_edge = 2
            subset_idx_cell = 5
            resize!(grid_ggdn.grid_subset, subset_idx_cell)  # Nodes, faces, x-aligned faces, y-aligned faces, cells
            grid_ggdn.grid_subset[subset_idx_node].dimension = 1
            grid_ggdn.grid_subset[subset_idx_node].identifier.name = "nodes"
            grid_ggdn.grid_subset[subset_idx_node].identifier.index = subset_idx_node
            grid_ggdn.grid_subset[subset_idx_node].identifier.description = "all points in the domain"
            grid_ggdn.grid_subset[subset_idx_edge].dimension = 2
            grid_ggdn.grid_subset[subset_idx_edge].identifier.name = "faces"
            grid_ggdn.grid_subset[subset_idx_edge].identifier.index = subset_idx_edge
            grid_ggdn.grid_subset[subset_idx_edge].identifier.description = "All edges in the domain"
            grid_ggdn.grid_subset[subset_idx_cell].dimension = 3
            grid_ggdn.grid_subset[subset_idx_cell].identifier.name = "cells"
            grid_ggdn.grid_subset[subset_idx_cell].identifier.index = subset_idx_cell
            grid_ggdn.grid_subset[subset_idx_cell].identifier.description = "all 2d cells in the domain"

            # Resizing objects to hold cell data
            # Should be fewer than this many points, but this way we won't under-fill
            resize!(o0.object, ncell * 4)  # Points
            resize!(o1.object, ncell * 4)  # Edges
            resize!(o2.object, ncell)  # Faces
            resize!(o3.object, ncell)  # Volumes

            for i = 1:(ny*nx*4)
                o0.object[i].geometry = [0.0, 0.0]
            end
            for i = 1:(ny*nx)
                o2.object[i].nodes = [0, 0, 0, 0]
            end

            resize!(dd.edge_profiles.ggd,nt)
            ggd = dd.edge_profiles.ggd
            for it in 1:nt
                if time_data_present
                    ggd[it].time = Float64.(times[it])
                end
                grp = getproperty(ggd[it], Symbol(group))
                qty = getproperty(grp, Symbol(quantity))
                resize!(qty, subset_idx_cell)
                qty[subset_idx_cell].grid_index = 1
                qty[subset_idx_cell].grid_subset_index = subset_idx_cell
                #resize!(qty[subset_idx_cell].values, nx*ny)
                qty[subset_idx_cell].values = zeros(Float64, nx*ny)
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
                        i_existing = search_points(dd, crx[icorner, iy, ix], cry[icorner, iy, ix])[1]
                        if i_existing == 0
                            o0.object[j].geometry = [crx[icorner, iy, ix], cry[icorner, iy, ix]]
                            o2.object[ic].nodes[icorner] = j
                            j += 1
                        else
                            o2.object[ic].nodes[icorner] = i_existing[1]
                        end
                    end
                    # Adding the data quantity of each cell in ggd.group.quantity[5].values[:]
                    for it=1:nt
                        grp = getproperty(ggd[it], Symbol(group))
                        qty = getproperty(grp, Symbol(quantity))
                        if time_data_present
                            qty[subset_idx_cell].values[ic] = values[it, iy, ix]  # one per cell
                        else
                            qty[subset_idx_cell].values[ic] = values[iy, ix]  # one per cell
                        end
                    end
                end  # for ix
            end  # for iy
        end  # for space_number
    end  # for grid_number
    return dd
end

end # module SOLPS2IMAS
