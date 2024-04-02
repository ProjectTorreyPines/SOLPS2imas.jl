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

function search_points(nodes, r, z; tol=0)
    n = length(r)
    indices = zeros(Int, n)
    # If an index remains at 0, it means the point in question was not found
    for j ∈ 1:n
        for i ∈ eachindex(nodes)
            rn, zn = getfield(nodes[i], :geometry)
            # zn = nodes[i].geometry[2]
            if abs(rn - r[j]) <= tol && abs(zn - z[j]) <= tol
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
        edges_ii_nodes = getfield(edges[ii], :nodes)
        if edge_nodes[1] == edges_ii_nodes[1] && edge_nodes[2] == edges_ii_nodes[2]
            return ii
        elseif edge_nodes[2] == edges_ii_nodes[1] &&
               edge_nodes[1] == edges_ii_nodes[2]
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

"""
    neighbour_inds(ic; nx, ny, leftcut, rightcut, topcut, bottomcut)

Returns indices of neighbours of cell with linear index ic. This function uses the SOLPS
grid generation algorithm to determine the neighbours. However, SOLPS geometry file
actually provides the neighbor indices directly. Thus, this function is not used in the
code anywhere but is kept here for reference.
"""
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

"""
    get_neighbour_inds(ic, gmtry, it)

Returns indices of neighbours of cell with linear index ic. This function uses the SOLPS
geometry file to determine the neighbours by using matrices named as leftix, rightix,
topix, bottomix, leftiy, rightiy, topiy, and bottomiy.
"""
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

"""
    attach_neightbours(cells, edges, gmtry, it)

This function attaches neighbours to each boundary of each cell and each boundary of
each edge. This function uses the SOLPS geometry file to determine the neighbours.
"""
function attach_neightbours(cells, edges, gmtry, it)
    for (ic, cell) ∈ enumerate(cells)
        for neighbour_ind ∈ get_neighbour_inds(ic, gmtry, it)
            for boundary ∈ cell.boundary
                for neighbour_boundary ∈ cells[neighbour_ind].boundary
                    if getfield(boundary, :index) ==
                       getfield(neighbour_boundary, :index) &&
                       neighbour_ind ∉ getfield(boundary, :neighbours)
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
                for edge_bnd ∈ getfield(edges[edge_ind], :boundary)
                    for neighbour_edge_bnd ∈ edges[neighbour_edge_ind].boundary
                        if getfield(edge_bnd, :index) ==
                           getfield(neighbour_edge_bnd, :index) &&
                           neighbour_edge_ind ∉ getfield(edge_bnd, :neighbours)
                            append!(edge_bnd.neighbours, neighbour_edge_ind)
                        end
                    end
                end
            end
        end
    end
end
