"""
in_core(; ix, iy, topcut, bottomcut, leftcut, rightcut)::Bool

Returns true if cell indexed `ix`, `iy` lie inside the core.
"""
function in_core(; ix, iy, topcut, bottomcut, leftcut, rightcut)::Bool
    return bottomcut + 1 < iy < topcut + 2 && leftcut + 1 < ix < rightcut + 2
end

"""
    in_sol(; iy, topcut, kwargs...)::Bool

Returns true if cell indexed `ix`, `iy` lie inside the SOL.
"""
in_sol(; iy, topcut, kwargs...)::Bool = topcut + 1 < iy

"""
    in_idr(; ix, iy, topcut, bottomcut, leftcut, kwargs...)::Bool

Returns true if cell indexed `ix`, `iy` lie inside the inner divertor region.
"""
function in_idr(; ix, iy, topcut, bottomcut, leftcut, kwargs...)::Bool
    return bottomcut + 1 < iy < topcut + 2 && ix < leftcut + 2
end

"""
    in_odr(; ix, iy, topcut, bottomcut, rightcut, kwargs...)::Bool

Returns true if cell indexed `ix`, `iy` lie inside the outer divertor region.
"""
function in_odr(; ix, iy, topcut, bottomcut, rightcut, kwargs...)::Bool
    return bottomcut + 1 < iy < topcut + 2 && rightcut + 1 < ix
end

"""
    is_x_aligned(;boundary_ind)::Bool

x\\_aligned edges will have odd `boundary_ind` based on chosen order of numbering them.
"""
is_x_aligned(; boundary_ind)::Bool = mod(boundary_ind, 2) == 1

"""
    is_y_aligned(; boundary_ind)::Bool

y\\_aligned edges will have even `boundary_ind` based on chosen order of numbering them.
"""
is_y_aligned(; boundary_ind)::Bool = mod(boundary_ind, 2) == 0

"""
    is_core_cut(;
        ix,
        iy,
        cells,
        nx,
        boundary_ind,
        topcut,
        bottomcut,
        leftcut,
        rightcut,
    )::Bool

Returns true if `boundary_ind` of a cell at `ix`, `iy` is on core\\_cut (Y-aliged edge).
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
)::Bool
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
    is_pfr_cut(;
        ix,
        iy,
        cells,
        nx,
        boundary_ind,
        topcut,
        bottomcut,
        leftcut,
        rightcut,
    )::Bool

Returns true if `boundary_ind` of a cell at `ix`, `iy` is on core_cut (y-aliged edge).
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
)::Bool
    if bottomcut + 1 < iy < topcut + 2 && ix == leftcut + 1 && mod(boundary_ind, 2) == 0
        ixr = rightcut + 2
        this_cell = cells[xytoc(ix, iy; nx=nx)]
        other_side_cell_ind = xytoc(ixr, iy; nx=nx)  # Cell on the other side of pfr cut
        return other_side_cell_ind ∈ this_cell.boundary[boundary_ind].neighbours
    end
    return false
end

"""
    is_outer_throat(; ix, iy, boundary_ind, topcut, rightcut, kwargs...)::Bool

Returns true if `boundary_ind` of a cell at `ix`, `iy` is on outer throat.
"""
function is_outer_throat(; ix, iy, boundary_ind, topcut, rightcut, kwargs...)::Bool
    return topcut + 1 < iy && ix == rightcut + 1 && boundary_ind == 2
end

"""
    is_inner_throat(; ix, iy, boundary_ind, topcut, leftcut, kwargs...)::Bool

Returns true if `boundary_ind` of a cell at `ix`, `iy` is on outer throat.
"""
function is_inner_throat(; ix, iy, boundary_ind, topcut, leftcut, kwargs...)::Bool
    return topcut + 1 < iy && ix == leftcut + 2 && boundary_ind == 4
end

"""
    is_outer_midplane(; ix, jxa, boundary_ind)

Returns true if `boundary_ind` of a cell at `ix`, `iy` is on outer midplane.
"""
function is_outer_midplane(; ix, iy, jxa, boundary_ind, topcut, kwargs...)::Bool
    # Note: USING CONVENTION to mark bottom edge of the midplane cell as midplane
    return ix == jxa && boundary_ind == 2
end

"""
    is_inner_midplane(; ix, jxa, boundary_ind)

Returns true if `boundary_ind` of a cell at `ix`, `iy` is on outer midplane.
"""
function is_inner_midplane(; ix, iy, jxi, boundary_ind, topcut, kwargs...)::Bool
    # Note: USING CONVENTION to mark bottom edge of the midplane cell as midplane
    return ix == jxi && boundary_ind == 4
end

"""
    is_outer_target(; ix, nx, boundary_ind)

Returns true if `boundary_ind` of a cell at `ix`, `iy` is on outer target.
"""
is_outer_target(; ix, nx, boundary_ind) = ix == nx && boundary_ind == 2

"""
    is_inner_target(; ix, boundary_ind, kwargs...)::Bool

Returns true if `boundary_ind` of a cell at `ix`, `iy` is on inner target.
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

Returns true if `boundary_ind` of a cell at `ix`, `iy` is on core boundary (central
blank spot boundary).
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
    is_separatrix(; iy, boundary_ind, topcut, kwargs...)::Bool

Returns true if `boundary_ind` of a cell at `ix`, `iy` is on separatrix.
"""
function is_separatrix(; iy, boundary_ind, topcut, kwargs...)::Bool
    return topcut + 2 == iy && boundary_ind == 1
end

"""
    get_xpoint_nodes(
        gmtry::Dict{String, Dict{String, Any}},
    )::Vector{Vector{Vector{Float64}}}

Limited to finding first x-point for now. Returns x-point (r, z) for each time index of
grid\\_ggd for the first x-point only. Thus second index correspond to the rank of x-point
which is always 1 from output of this function for now.
"""
function get_xpoint_nodes(
    gmtry::Dict{String, Dict{String, Any}},
)::Vector{Vector{Vector{Float64}}}
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
    candidate_nodes = [
        [
            [
                [crx[it, icorner, iy, ix], cry[it, icorner, iy, ix]] for icorner ∈ 1:4
            ] for (iy, ix) ∈ xpcells
        ] for it ∈ 1:nt
    ]
    # Find the node that is common among all the cells
    xpoint_nodes = [intersect(candidate_nodes[it]...) for it ∈ 1:nt]
    return xpoint_nodes
end
