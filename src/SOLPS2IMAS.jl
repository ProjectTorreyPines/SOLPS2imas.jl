module SOLPS2IMAS

# using DelimitedFiles
using IMASDD
const OMAS = IMASDD

export try_omas
export populate_grid_ggd
export generate_test_data
export read_b2_output

function try_omas()
    println("it's the omas function")
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

function read_b2_output(filename)
    # x = readdlm(filename, ' ', Float64, '\n', header=true)  # Doesn't handle the text lines at the start of each array
    contents = Dict()
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
            # Cleanup previous array if applicable
            if tag != ""
                if arraysize == nx * ny
                    # println(tag, " ", arraysize, " ", nx, " ", ny, " ", size(contents[tag]))
                    # print(contents[tag])
                    contents[tag] = reshape(contents[tag], (ny, nx))
                elseif arraysize == nx * ny * ns
                    # If ns == 2, then r,z vector arrays can't be distinguished
                    # from species-dependent quantities by their shapes. But
                    # they get treated the same way, so it's okay.
                    contents[tag] = reshape(contents[tag], (ns, ny, nx))
                elseif arraysize == nx * ny * 2
                    contents[tag] = reshape(contents[tag], (2, ny, nx))
                elseif arraysize == nx * ny * 2 * ns
                    contents[tag] = reshape(contents[tag], (ns, 2, ny, nx))
                end
            end
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

            if tag == "nx,ny,ns"
                nx, ny, ns = array_line
                nx += 2  # Account for guard cells
                ny += 2  # Account for guard cells
            end
        end
    end
    return contents
end

function populate_grid_ggd(nx::Int64, ny::Int64)
    ncell = nx * ny
    dd = OMAS.dd()
    println("another fun function!!!!!11!!!!!")
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

    resize!(o0.object, ncell * 4)  # Points
    resize!(o1.object, ncell * 4)  # Edges
    resize!(o2.object, ncell)  # Faces
    resize!(o3.object, ncell)  # Volumes
    for iy = 1:ny
        for ix = 1:nx
            ic::Int = (iy - 1) * nx + ix
            resize!(o2.object[ic].boundary, 4)
        end # for ix
    end # for iy

    return nothing
end

end # module SOLPS2IMAS
