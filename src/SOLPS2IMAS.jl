module SOLPS2IMAS

using Revise
using OMAS: OMAS
using NCDatasets: Dataset, dimnames
using YAML: load_file as YAML_load_file
import GGDUtils: add_subset_element!, get_grid_subset_with_index, get_subset_boundary,
    get_subset_space, subset_do
export try_omas
export generate_test_data
export read_b2_output
export search_points
export solps2imas

function try_omas()
    ids = OMAS.dd()
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
    θₓ = 8 * π / 6
    θᵢ = θₓ - π / 4
    θₒ = θₓ + π / 4
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
    for ix::Int ∈ 1:nx
        for iy::Int ∈ 1:ny
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
                cell_centers_r[iy, ix] =
                    xpoint_r + cos(θᵢ) * dx_il * (lcut - ix) + sin(θᵢ) * dy * (iy - sep)
                cell_centers_z[iy, ix] =
                    xpoint_z + sin(θᵢ) * dx_il * (lcut - ix) + cos(θᵢ) * dy * (iy - sep)
            else
                cell_centers_r[iy, ix] =
                    xpoint_r + cos(θₒ) * dx_ol * (ix - rcut) + sin(θₒ) * dy * (iy - sep)
                cell_centers_z[iy, ix] =
                    xpoint_z + sin(θₒ) * dx_ol * (ix - rcut) + cos(θₒ) * dy * (iy - sep)
            end
        end
    end
    # Sanitize
    rmin = 0.0
    rmax = 3.0
    zmin = -3.0
    zmax = 3.0
    for i::Int ∈ 1:nx*ny
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
    for key ∈ keys(ds.dim)
        ret_dict["dim"][key] = ds.dim[key]
    end
    for key ∈ keys(ds)
        if key != "ntstep"
            d = dimnames(ds[key])
            permute = [
                y for
                y ∈ [findfirst(x -> x == dimord, d) for dimord ∈ dim_order] if
                y !== nothing
            ]
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

function read_b2mn_output(filename)
    lines = open(filename) do f
        return readlines(f)
    end
    contents = Dict()
    for line ∈ lines
        if startswith(line, "'")
            # Ignore comments and remove spaces
            line = strip(split(line, "#")[1], [' '])
            splits = split(line, "'"; keepempty=false)
            contents[splits[1]] = parse(Float64, splits[end])
        end
    end
    return contents
end

function read_b2_output(filename)
    if cmp(splitext(filename)[2], ".nc") == 0
        return read_b2time_output(filename)
    end

    contents = Dict()
    array_sizes = Dict()
    ret_dict = Dict()
    lines = open(filename) do f
        return readlines(f)
    end
    nx = 0
    ny = 0
    ns = 0
    tag = ""
    arraysize = 0
    arraytype = nothing
    j = 1
    for l ∈ lines
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
                array_line = [parse(Int, ss) for ss ∈ split(l)]
                array_inc = size(array_line)[1]
            elseif arraytype == "real"
                array_line = [parse(Float64, ss) for ss ∈ split(l)]
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
        error(
            "nx,ny (b2fgmtry) or nx,ny,ns (b2fstate) must be present in b2 output file",
        )
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
    for k ∈ keys(gmtry)
        # The 4 fields of bb are poloidal, radial, toroidal, and total magnetic field
        # according to page 212 of D. Coster, "SOLPS-ITER [manual]" (2019)
        # The 4 fields in crx and cry are the corners of each grid cell.
        if k ∈ ["crx", "cry", "bb"]
            ret_dict["data"][k] =
                permutedims(reshape(gmtry[k], (nx, ny, 4, 1)), (4, 3, 2, 1))
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
    ret_dict["dim"]["nx_no_guard"],
    ret_dict["dim"]["ny_no_guard"],
    ret_dict["dim"]["ns"] = state["nx,ny,ns"]
    # includes guard cells
    nx = ret_dict["dim"]["nx"] = ret_dict["dim"]["nx_no_guard"] + 2
    ny = ret_dict["dim"]["ny"] = ret_dict["dim"]["ny_no_guard"] + 2
    ns = ret_dict["dim"]["ns"]
    ndir = ret_dict["dim"]["ndir"] = 2
    # Adding placeholder timestamp
    ret_dict["dim"]["time"] = 1
    ret_dict["data"]["timesa"] = [0.0]
    for k ∈ keys(state)
        l = length(state[k])
        if l == nx * ny
            ret_dict["data"][k] = permutedims(reshape(state[k], (nx, ny, 1)), (3, 2, 1))
        elseif l == nx * ny * ns
            ret_dict["data"][k] =
                permutedims(reshape(state[k], (nx, ny, ns, 1)), (4, 3, 2, 1))
        elseif l == nx * ny * ndir
            ret_dict["data"][k] =
                permutedims(reshape(state[k], (nx, ny, ndir, 1)), (4, 3, 2, 1))
        elseif l == nx * ny * ndir * ns
            ret_dict["data"][k] =
                permutedims(reshape(state[k], (nx, ny, ndir, ns, 1)), (5, 4, 3, 2, 1))
        elseif l == ns
            ret_dict["data"][k] = permutedims(reshape(state[k], (ns, 1)), (2, 1))
        elseif k ∉ keys(ret_dict["dim"])
            ret_dict["data"][k] = state[k]
        end
    end
    return ret_dict
end

"""
    xytoc(ix, iy; nx)

Converts SOLPS indices for crx, cry (ix, iy) that go from 1:nx, 1:ny
into the linear index ic used in IMAS for corresponding cells
"""
function xytoc(ix, iy; nx)
    ic::Int = (iy - 1) * nx + ix
    return ic
end

"""
    ctoxy(ic; nx)

Inverse of xytoc
"""
function ctoxy(ic; nx)
    ix::Int = mod(ic - 1, nx) + 1
    iy::Int = (ic - 1) ÷ nx + 1
    return ix, iy
end

"""
    data_xytoc(data; nx)

Flattens 2d data given on cell indices ix, iy into 1d data on linear index ic. ic is
calculated using xytoc function. Data is assumed to have dimensions (ny, nx)
"""
function data_xytoc(data; nx)
    flat_data = Array{eltype(data)}(undef, length(data))
    for ix ∈ axes(data, 2)
        for iy ∈ axes(data, 1)
            flat_data[xytoc(ix, iy; nx=nx)] = data[iy, ix]
        end
    end
    return flat_data
end

"""
in_core(; ix, iy, topcut, bottomcut, leftcut, rightcut)

Returns true if cell indexed ix, iy lie inside the core
"""
function in_core(; ix, iy, topcut, bottomcut, leftcut, rightcut)
    return bottomcut + 1 < iy < topcut + 2 && leftcut + 1 < ix < rightcut + 2
end

"""
    in_sol(; iy, topcut, kwargs...)

Returns true if cell indexed ix, iy lie inside the SOL
"""
in_sol(; iy, topcut, kwargs...) = topcut + 1 < iy

"""
    in_idr(; ix, iy, topcut, bottomcut, leftcut, kwargs...)

Returns true if cell indexed ix, iy lie inside the inner divertor region
"""
function in_idr(; ix, iy, topcut, bottomcut, leftcut, kwargs...)
    return bottomcut + 1 < iy < topcut + 2 && ix < leftcut + 2
end

"""
    in_odr(; ix, iy, topcut, bottomcut, rightcut, kwargs...)

Returns true if cell indexed ix, iy lie inside the outer divertor region
"""
function in_odr(; ix, iy, topcut, bottomcut, rightcut, kwargs...)
    return bottomcut + 1 < iy < topcut + 2 && rightcut + 1 < ix
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

"""
    is_x_aligned(;boundary_ind)

y_aligned edges will have odd boundary_ind based on chosen order of numbering them
"""
is_x_aligned(; boundary_ind) = mod(boundary_ind, 2) == 1

"""
    is_y_aligned(; boundary_ind)

y_aligned edges will have even boundary_ind based on chosen order of numbering them
"""
is_y_aligned(; boundary_ind) = mod(boundary_ind, 2) == 0

"""
is_core_cut(; ix, iy, nx, cells, boundary_ind, topcut, bottomcut, leftcut, rightcut)

Returns true if boundary_ind of a cell at ix, iy is on core_cut (Y-aliged edge)
"""
function is_core_cut(;
    ix,
    iy,
    cells,
    nx,
    boundary_ind,
    topcut,
    bottomcut,
    leftcut,
    rightcut,
)
    if bottomcut + 1 < iy < topcut + 2 && ix == leftcut + 2 && mod(boundary_ind, 2) == 0
        ixr = rightcut + 1
        this_cell = cells[xytoc(ix, iy; nx=nx)]
        # Cell on the other side of core cut
        other_side_cell_ind = xytoc(ixr, iy; nx=nx)
        return other_side_cell_ind ∈ this_cell.boundary[boundary_ind].neighbours
    end
    return false
end

"""
is_pfr_cut(; ix, iy, nx, cells, boundary_ind, topcut, bottomcut, leftcut, rightcut)

Returns true if boundary_ind of a cell at ix, iy is on core_cut (Y-aliged edge)
"""
function is_pfr_cut(;
    ix,
    iy,
    cells,
    nx,
    boundary_ind,
    topcut,
    bottomcut,
    leftcut,
    rightcut,
)
    if bottomcut + 1 < iy < topcut + 2 && ix == leftcut + 1 && mod(boundary_ind, 2) == 0
        ixr = rightcut + 2
        this_cell = cells[xytoc(ix, iy; nx=nx)]
        other_side_cell_ind = xytoc(ixr, iy; nx=nx)  # Cell on the other side of pfr cut
        return other_side_cell_ind ∈ this_cell.boundary[boundary_ind].neighbours
    end
    return false
end

"""
    is_outer_throat(; ix, iy, boundary_ind, topcut, rightcut, kwargs...)

Returns true if boundary_ind of a cell at ix, iy is on outer throat
"""
function is_outer_throat(; ix, iy, boundary_ind, topcut, rightcut, kwargs...)
    return topcut + 1 < iy && ix == rightcut + 1 && boundary_ind == 2
end

"""
    is_inner_throat(; ix, iy, boundary_ind, topcut, leftcut, kwargs...)

Returns true if boundary_ind of a cell at ix, iy is on outer throat
"""
function is_inner_throat(; ix, iy, boundary_ind, topcut, leftcut, kwargs...)
    return topcut + 1 < iy && ix == leftcut + 2 && boundary_ind == 4
end

"""
    is_outer_midplane(; ix, jxa, boundary_ind)

Returns true if boundary_ind of a cell at ix, iy is on outer midplane
"""
function is_outer_midplane(; ix, iy, jxa, boundary_ind, topcut, kwargs...)
    # Note: USING CONVENTION to mark bottom edge of the midplane cell as midplane
    return ix == jxa && boundary_ind == 2
end

"""
    is_inner_midplane(; ix, jxa, boundary_ind)

Returns true if boundary_ind of a cell at ix, iy is on outer midplane
"""
function is_inner_midplane(; ix, iy, jxi, boundary_ind, topcut, kwargs...)
    # Note: USING CONVENTION to mark bottom edge of the midplane cell as midplane
    return ix == jxi && boundary_ind == 4
end

"""
    is_outer_target(; ix, nx, boundary_ind)

Returns true if boundary_ind of a cell at ix, iy is on outer target
"""
is_outer_target(; ix, nx, boundary_ind) = ix == nx && boundary_ind == 2

"""
    is_inner_target(; ix, boundary_ind, kwargs...)

Returns true if boundary_ind of a cell at ix, iy is on inner target
"""
is_inner_target(; ix, boundary_ind) = ix == 1 && boundary_ind == 4

"""
    is_core_boundary(;
    ix,
    iy,
    boundary_ind,
    bottomcut,
    leftcut,
    rightcut,
    kwargs...,

)

Returns true if boundary_ind of a cell at ix, iy is on core boundary (central blank
spot boundary)
"""
function is_core_boundary(;
    ix,
    iy,
    boundary_ind,
    bottomcut,
    leftcut,
    rightcut,
    kwargs...,
)
    return bottomcut + 2 == iy && leftcut + 1 < ix < rightcut + 2 && boundary_ind == 1
end

"""
    is_separatix(; iy, boundary_ind, topcut, kwargs...)

Returns true if boundary_ind of a cell at ix, iy is on separatix
"""
function is_separatix(; iy, boundary_ind, topcut, kwargs...)
    return topcut + 2 == iy && boundary_ind == 1
end

"""
    add_subset_element!(
    subset,
    sn,
    dim,
    index::Int,
    in_subset=(x...) -> true;
    kwargs...,

)

Adds the geometric element in subset object (assumed to be resized already) at element
dd index dd_ind, with space number sn, dimension dim, index. To determine,
if the element should be added or not, a function in_subset can be provided that gets
the arguments (kwargs...). These functions will be in_core, in_sol etc as difined above.
"""
function add_subset_element!(
    subset,
    sn,
    dim,
    index::Int,
    in_subset=(x...) -> true;
    kwargs...,
)
    if in_subset(; kwargs...)
        dd_ind = length(subset.element) + 1
        resize!(subset.element, dd_ind)
        resize!(subset.element[dd_ind].object, 1)
        subset.element[dd_ind].object[1].space = sn
        subset.element[dd_ind].object[1].dimension = dim
        subset.element[dd_ind].object[1].index = index
    end
end

"""
    add_subset_element!(
    subset,
    sn,
    dim,
    index::Vector{Int},
    in_subset=(x...) -> true;
    kwargs...,

)

Overloaded to work differently (faster) with list of indices to be added.
"""
function add_subset_element!(
    subset,
    sn,
    dim,
    index::Vector{Int},
    in_subset=(x...) -> true;
    kwargs...,
)
    if in_subset(; kwargs...)
        dd_start_ind = length(subset.element) + 1
        resize!(subset.element, length(subset.element) + length(index))
        dd_stop_ind = length(subset.element)
        for (ii, dd_ind) ∈ enumerate(dd_start_ind:dd_stop_ind)
            resize!(subset.element[dd_ind].object, 1)
            subset.element[dd_ind].object[1].space = sn
            subset.element[dd_ind].object[1].dimension = dim
            subset.element[dd_ind].object[1].index = index[ii]
        end
    end
end

"""
    get_xpoint_nodes(gmtry)

Limited to finding first x-point for now.
"""
function get_xpoint_nodes(gmtry)
    crx = gmtry["data"]["crx"]
    cry = gmtry["data"]["cry"]
    nt = gmtry["dim"]["time"]
    # Find cells around x-point
    xpcells = [(gmtry["data"]["topcut"][1] + 1, gmtry["data"]["leftcut"][1] + 1),
        (gmtry["data"]["topcut"][1] + 1, gmtry["data"]["leftcut"][1] + 2),
        (gmtry["data"]["topcut"][1] + 2, gmtry["data"]["leftcut"][1] + 1),
        (gmtry["data"]["topcut"][1] + 2, gmtry["data"]["leftcut"][1] + 2),
        (gmtry["data"]["topcut"][1] + 1, gmtry["data"]["rightcut"][1] + 1),
        (gmtry["data"]["topcut"][1] + 1, gmtry["data"]["rightcut"][1] + 2),
        (gmtry["data"]["topcut"][1] + 2, gmtry["data"]["rightcut"][1] + 1),
        (gmtry["data"]["topcut"][1] + 2, gmtry["data"]["rightcut"][1] + 2)]
    # Get list of all nodes in these cells
    candidate_nodes = []
    resize!(candidate_nodes, nt)
    for it ∈ 1:nt
        candidate_nodes[it] = [
            [
                [crx[it, icorner, iy, ix], cry[it, icorner, iy, ix]] for icorner ∈ 1:4
            ] for (iy, ix) ∈ xpcells
        ]
    end
    xpoint_nodes = []
    resize!(xpoint_nodes, nt)
    # Find the node that is common among all the cells
    for it ∈ 1:nt
        xpoint_nodes[it] = intersect(candidate_nodes[it]...)
    end
    return xpoint_nodes
end

function search_points(nodes, r, z)
    n = length(r)
    indices = zeros(Int, n)
    # If an index remains at 0, it means the point in question was not found
    for j ∈ 1:n
        for i ∈ eachindex(nodes)
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
    search_edges(edges, edge_nodes)

search if an edge with nodes as edge_nodes already exists
"""
function search_edges(edges, edge_nodes)
    for ii ∈ eachindex(edges)
        if edge_nodes[2] == edges[ii].nodes[1] && edge_nodes[1] == edges[ii].nodes[2]
            return ii
        elseif edge_nodes[2] == edges[ii].nodes[1] &&
               edge_nodes[1] == edges[ii].nodes[2]
            return ii
        end
    end
    return 0
end

"""
    distance_between_nodes(nodes, node_inds)

Return distance between two node indices
"""
function distance_between_nodes(nodes, node_inds)
    return √(sum((nodes[node_inds[1]].geometry - nodes[node_inds[2]].geometry) .^ 2))
end

function neighbour_inds(ic; nx, ny, leftcut, rightcut, topcut, bottomcut)
    ix, iy = ctoxy(ic; nx=nx)
    neighbour_x_inds = []
    neighbour_y_inds = []
    if ix > 1
        if bottomcut + 1 < iy ≤ topcut + 1
            if ix == rightcut + 2  # left most outter divertor region
                append!(neighbour_x_inds, leftcut + 1)
            elseif ix == leftcut + 2  # left most core region
                append!(neighbour_x_inds, rightcut + 1)
            else
                append!(neighbour_x_inds, ix - 1)
            end
        else
            append!(neighbour_x_inds, ix - 1)
        end
    end
    if ix < nx
        if bottomcut + 1 < iy ≤ topcut + 1
            if ix == leftcut + 1  # right most inner divertor regio
                append!(neighbour_x_inds, rightcut + 2)
            elseif ix == rightcut + 1  # right most core region
                append!(neighbour_x_inds, leftcut + 2)
            else
                append!(neighbour_x_inds, ix + 1)
            end
        else
            append!(neighbour_x_inds, ix + 1)
        end
    end
    if iy > 1
        append!(neighbour_y_inds, iy - 1)
    end
    if iy < ny
        append!(neighbour_y_inds, iy + 1)
    end

    neighbour_inds = []
    for x_ind ∈ neighbour_x_inds
        append!(neighbour_inds, xytoc(x_ind, iy; nx=nx))
    end
    for y_ind ∈ neighbour_y_inds
        append!(neighbour_inds, xytoc(ix, y_ind; nx=nx))
    end
    return neighbour_inds
end

function get_neighbour_inds(ic, gmtry, it)
    nx = gmtry["dim"]["nx"]
    ny = gmtry["dim"]["ny"]
    ix, iy = ctoxy(ic; nx=nx)
    neighbour_inds = []
    # println(ix, ", ", iy)
    for neighbour ∈ ["left", "right", "top", "bottom"]
        nix = gmtry["data"][neighbour*"ix"][it, iy, ix] + 2
        niy = gmtry["data"][neighbour*"iy"][it, iy, ix] + 2
        # println(neighbour, ": ", nix, ", ", niy)
        if 1 ≤ nix ≤ nx && 1 ≤ niy ≤ ny
            append!(neighbour_inds, xytoc(nix, niy; nx=nx))
        end
    end
    return neighbour_inds
end

function attach_neightbours(cells, edges, gmtry, it)
    for (ic, cell) ∈ enumerate(cells)
        for neighbour_ind ∈ get_neighbour_inds(ic, gmtry, it)
            for boundary ∈ cell.boundary
                for neighbour_boundary ∈ cells[neighbour_ind].boundary
                    if boundary.index == neighbour_boundary.index &&
                       neighbour_ind ∉ boundary.neighbours
                        append!(boundary.neighbours, neighbour_ind)
                    end
                end
            end
        end
    end
    for (ic, cell) ∈ enumerate(cells)
        for edge_ind ∈ [bnd.index for bnd ∈ cell.boundary]
            neighbour_edge_inds = [bnd.index for bnd ∈ cell.boundary]
            for neighbour_ind ∈ get_neighbour_inds(ic, gmtry, it)
                union!(
                    neighbour_edge_inds,
                    [bnd.index for bnd ∈ cells[neighbour_ind].boundary],
                )
            end
            setdiff!(neighbour_edge_inds, edge_ind)
            for neighbour_edge_ind ∈ neighbour_edge_inds
                for edge_bnd ∈ edges[edge_ind].boundary
                    for neighbour_edge_bnd ∈ edges[neighbour_edge_ind].boundary
                        if edge_bnd.index == neighbour_edge_bnd.index &&
                           neighbour_edge_ind ∉ edge_bnd.neighbours
                            append!(edge_bnd.neighbours, neighbour_edge_ind)
                        end
                    end
                end
            end
        end
    end
end

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
YAML file. Returns data in OMAS.dd datastructure.
"""
function solps2imas(b2gmtry, b2output, gsdesc, b2mn=nothing; load_bb=false)
    # Initialize an empty OMAS data structre
    ids = OMAS.dd()

    # Setup the grid first
    gmtry = read_b2_output(b2gmtry)

    jxi = jxa = nothing
    if !isnothing(b2mn)
        mn = read_b2mn_output(b2mn)
        if "b2mwti_jxa" ∈ keys(mn)
            jxa = Int(mn["b2mwti_jxa"])
        end
        if "b2mwti_jxi" ∈ keys(mn)
            jxi = Int(mn["b2mwti_jxa"])
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
