module SOLPS2IMAS

using IMASDD
const OMAS = IMASDD

export try_omas
export populate_grid_ggd

function try_omas()
    println("it's the omas function")
    dd = OMAS.dd()
    resize!(dd.equilibrium.time_slice, 1)
    dd.equilibrium.time_slice[1].profiles_1d.psi = [0.0, 1.0, 2.0, 3.9]
    return nothing
end

function populate_grid_ggd()
    dd = OMAS.dd()
    resize!(dd.edge_profiles, 1)
    resize!(dd.edge_profiles.grid_ggd, 1)
    dd.edge_profiles.grid_ggd[1].identifier.name = "Sven"
    return nothing
end

end # module SOLPS2IMAS
