"""
    xytoc(ix::Int, iy::Int; nx::Int)::Int

Converts SOLPS indices for crx, cry (`ix`, `iy`) that go from 1:`nx`, 1:`ny`
into the linear index `ic` used in IMAS for corresponding cells.
"""
function xytoc(ix::Int, iy::Int; nx::Int)::Int
    ic::Int = (iy - 1) * nx + ix
    return ic
end

"""
    ctoxy(ic::Int; nx::Int)::Tuple{Int, Int}

Inverse of xytoc.
"""
function ctoxy(ic::Int; nx::Int)::Tuple{Int, Int}
    ix::Int = mod(ic - 1, nx) + 1
    iy::Int = (ic - 1) ÷ nx + 1
    return ix, iy
end

"""
    data_xytoc(data::Matrix{Float64}; nx::Int)::Vector{Float64}

Flattens 2d data given on cell indices `ix`, `iy` into 1d data on linear index `ic`.
`ic` is calculated using xytoc function. Data is assumed to have dimensions (ny, `nx`)
where ny is not required in this conversion.
"""
function data_xytoc(data::Matrix{Float64}; nx::Int)::Vector{Float64}
    flat_data = Array{eltype(data)}(undef, length(data))
    for ix ∈ axes(data, 2)
        for iy ∈ axes(data, 1)
            flat_data[xytoc(ix, iy; nx=nx)] = data[iy, ix]
        end
    end
    return flat_data
end

edges_nodes_type = IMASdd.IDSvector{
    IMASdd.edge_profiles__grid_ggd___space___objects_per_dimension___object{T},
} where {T}

"""
    search_point(
        nodes::IMASdd.IDSvector{IMASdd.edge_profiles__grid_ggd___space___objects_per_dimension___object{T}},
        r::Real,
        z::Real;
        tol::Float64=0.0,
    )::Int where {T}

Search if a point (`r`, `z`) is present in the nodes array. Here `nodes` is generally
available in `ids.edge_profiles.grid_ggd[:].space[:].objects_per_dimension[1].object`

If the point is not found, the function returns 0.
"""
function search_point(
    nodes::edges_nodes_type,
    r::Real,
    z::Real;
    tol::Float64=0.0,
)::Int
    for ii ∈ eachindex(nodes)
        rn, zn = getfield(nodes[ii], :geometry)
        if abs(rn - r) <= tol && abs(zn - z) <= tol
            return ii
        end
    end
    return 0
end

"""
    search_edges(
        edges::IMASdd.IDSvector{IMASdd.edge_profiles__grid_ggd___space___objects_per_dimension___object{T}},
        edge_nodes::Array{Int, 1}
    )::Int where {T}

Search if an edge with nodes as `edge_nodes` already exists in the `edges` array.
`edges` is generally available in
`ids.edge_profiles.grid_ggd[:].space[:].objects_per_dimension[2].object`

If the edge is not found, the function returns 0.
"""
function search_edges(edges::edges_nodes_type, edge_nodes::Array{Int, 1})::Int
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
    distance_between_nodes(nodes::edges_nodes_type, node_inds::Array{Int, 1})

Return distance between two nodes with indices `node_inds` in `nodes` array. `nodes` is
generally available in
`ids.edge_profiles.grid_ggd[:].space[:].objects_per_dimension[1].object`.
"""
function distance_between_nodes(nodes::edges_nodes_type, node_inds::Array{Int, 1})
    return √(sum((nodes[node_inds[1]].geometry - nodes[node_inds[2]].geometry) .^ 2))
end

"""
    neighbour_inds(
        ic::Int;
        nx::Int,
        ny::Int,
        leftcut::Int,
        rightcut::Int,
        topcut::Int,
        bottomcut::Int,
    )::Vector{Int}

(deprecated function)

Returns indices of neighbours of cell with linear index `ic`. This function uses the
SOLPS grid generation algorithm to determine the neighbours. However, SOLPS geometry
file actually provides the neighbor indices directly. Thus, this function is not used
in the code anywhere but is kept here for reference.
"""
function neighbour_inds(
    ic::Int;
    nx::Int,
    ny::Int,
    leftcut::Int,
    rightcut::Int,
    topcut::Int,
    bottomcut::Int,
)::Vector{Int}
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

    neighbour_inds = Vector{Int}(undef, 0)
    for x_ind ∈ neighbour_x_inds
        append!(neighbour_inds, xytoc(x_ind, iy; nx=nx))
    end
    for y_ind ∈ neighbour_y_inds
        append!(neighbour_inds, xytoc(ix, y_ind; nx=nx))
    end
    return neighbour_inds
end

"""
    get_neighbour_inds(
        ic::Int,
        gmtry::Dict{String, Dict{String, Any}},
        it::Int,
    )::Vector{Int}

Returns indices of neighbours of cell with linear index `ic`. This function uses the
parsed SOLPS geometry file to determine the neighbours by using matrices named as
leftix, rightix, topix, bottomix, leftiy, rightiy, topiy, and bottomiy.
"""
function get_neighbour_inds(
    ic::Int,
    gmtry::Dict{String, Dict{String, Any}},
    it::Int,
)::Vector{Int}
    nx = gmtry["dim"]["nx"]
    ny = gmtry["dim"]["ny"]
    ix, iy = ctoxy(ic; nx=nx)
    neighbour_inds = Vector{Int}(undef, 0)
    for neighbour ∈ ["left", "right", "top", "bottom"]
        nix = gmtry["data"][neighbour*"ix"][it, iy, ix] + 2
        niy = gmtry["data"][neighbour*"iy"][it, iy, ix] + 2
        if 1 ≤ nix ≤ nx && 1 ≤ niy ≤ ny
            append!(neighbour_inds, xytoc(nix, niy; nx=nx))
        end
    end
    return neighbour_inds
end

"""
    attach_neightbours!(
        cells::IMASdd.IDSvector{IMASdd.edge_profiles__grid_ggd___space___objects_per_dimension___object{T}},
        edges::IMASdd.IDSvector{IMASdd.edge_profiles__grid_ggd___space___objects_per_dimension___object{T}},
        gmtry::Dict{String, Dict{String, Any}},
        it::Int,
    ) where {T}

This function attaches neighbours to each boundary of each cell and each boundary of
each edge using the parsed SOLPS geometry file.
"""
function attach_neightbours!(
    cells::edges_nodes_type,
    edges::edges_nodes_type,
    gmtry::Dict{String, Dict{String, Any}},
    it::Int,
)
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
