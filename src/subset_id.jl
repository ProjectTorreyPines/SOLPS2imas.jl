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
