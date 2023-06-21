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
    println("another fun function!!!!!!!!!!")
    resize!(dd.edge_profiles.grid_ggd, 1)

    id = dd.edge_profiles.grid_ggd[1].identifier
    id.name = "Sven"
    id.index = 1
    id.description = "this is a grid"

    resize!(dd.edge_profiles.grid_ggd[1].space, 1)
    space = dd.edge_profiles.grid_ggd[1].space[1]
    space.identifier.name = "sp4ce"
    space.identifier.index = 1
    space.identifier.description = "The final frontier"
    space.geometry_type.name = "standard"  # I doubt this is needed
    space.geometry_type.index = 0  # 0 for standard, 1 for fourier. This is the important field
    space.geometry_type.description = "trying to hold a b2/solps mesh here"  # I doubt this is needed

    resize!(space.objects_per_dimension, 4)
    o0 = space.objects_per_dimension[1]  # 0D objects
    o1 = space.objects_per_dimension[2]  # 1D objects
    o2 = space.objects_per_dimension[3]  # 2D objects
    o3 = space.objects_per_dimension[4]  # 3D objects


    return nothing
end

end # module SOLPS2IMAS
