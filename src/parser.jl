export read_b2_output, read_b2mn_output, read_b2time_output # , read_b2_boundary_parameters

"""
    read_b2time_output(filename::String)::Dict{String, Dict{String, Any}}

Read time dependent b2 output file and return a dictionary with structure:
Dict("dim" => Dict{String, Any}, "data" => Dict{String, Any})
where "dim" contains the dimensions of the data and "data" contains the data itself,
with keys corresponding to the field names.

Supported SOLPS files as input via filename:

  - b2time.nc
"""
function read_b2time_output(filename::String)::Dict{String, Dict{String, Any}}
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
            scale = 1.0
            if "scale" in keys(ds[key].attrib)
                scale = ds[key].attrib["scale"]
            end
            try
                ret_dict["data"][key] = permutedims(Array(ds[key]), permute) .* scale
            catch e
                println("Error in reading ", key)
                showerror(stdout, e)
                println("Continuing by ignoring this field")
            end
        end
    end
    return ret_dict
end

"""
    read_b2mn_output(filename::String)::Dict{String, Any}

Read b2mn output file and store the quantities in a dictionary.

Supported SOLPS files as input via filename:

  - b2mn.dat
"""
function read_b2mn_output(filename::String)::Dict{String, Any}
    # Get list of integer fields
    d = readdlm("$(@__DIR__)/b2mn_int_fields.txt")
    int_fields = d[:, 1]
    # Get a dictionary of default defined fields
    def_int_fields =
        Dict(d[ii, 1] => d[ii, 2] for ii ∈ range(1, size(d)[1]) if isa(d[ii, 2], Int))
    lines = open(filename) do f
        return readlines(f)
    end
    contents = Dict()
    found_endphy = false
    for line ∈ lines
        # Remove all whitespace characters (taken from Base.isspace definition)
        line = strip(line, ['\t', '\n', '\v', '\f', '\r', ' ', '\u85', '\ua0'])
        if !found_endphy
            # Ignore all lines until *endphy
            if startswith(line, "*endphy")
                found_endphy = true
            end
            continue
        else
            if startswith(line, "'") || startswith(line, "\"")
                # Ignore comments that can start with #, !, or *
                line = split(line, "#"; keepempty=false)[1]
                line = split(line, "!"; keepempty=false)[1]
                line = split(line, "*"; keepempty=false)[1]
                # Replace all spaces with nothing, double quotes with single quotes
                line = replace(line, Base.isspace => "", "\"" => "'")
                # Now split with single quotes and discard empty strings
                name_value = split(line, "'"; keepempty=false)
                if length(name_value) == 0 # Case where key is commented inside quotes
                    continue
                end
                # Get key and value in lowercase
                key = lowercase(name_value[1])
                value = lowercase(name_value[2])
                try
                    value = parse(Float64, value)
                catch
                    value = parse(String, value)
                end
                if key in int_fields
                    value = Int(value)
                end
                contents[key] = value
            end
        end
    end
    for key ∈ keys(def_int_fields)
        if key ∉ keys(contents)
            contents[key] = def_int_fields[key]
        end
    end
    return contents
end

"""
    read_b2_output(filename::String)::Dict{String, Dict{String, Any}}

Read final state b2 output file (b2fstate or b2time.nc) or b2fgmtry file and return a
dictionary with structure:
Dict("dim" => Dict{String, Any}, "data" => Dict{String, Any})
where "dim" contains the dimensions of the data and "data" contains the data itself,
with keys corresponding to the field names.

Supported SOLPS files as input via filename:

  - b2fstate
  - b2fstati
  - b2time.nc
  - b2fgmtry
"""
function read_b2_output(filename::String)::Dict{String, Dict{String, Any}}
    if cmp(splitext(filename)[2], ".nc") == 0
        return read_b2time_output(filename)
    end

    contents = Dict{String, Any}()
    array_sizes = Dict{String, Any}()
    lines = open(filename) do f
        return readlines(f)
    end
    tag = ""
    arraysize = 0
    arraytype = nothing
    j = 1
    for l ∈ lines
        if startswith(l, "*cf:")
            j = 1  # Reset intra-array element counter
            _, arraytype, arraysize, tag = split(l)
            tag = String(tag)
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

function extract_geometry(gmtry::Dict{String, Any})::Dict{String, Dict{String, Any}}
    ret_dict = Dict("dim" => Dict{String, Any}(), "data" => Dict{String, Any}())
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

function extract_state_quantities(
    state::Dict{String, Any},
)::Dict{String, Dict{String, Any}}
    ret_dict = Dict("dim" => Dict{String, Any}(), "data" => Dict{String, Any}())
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
    read_b2_boundary_parameters(filename::String)::Dict{String, Any}

Reads and interprets the b2.boundary.parameters file from the SOLPS input deck.
This file has boundary conditions like power crossing into the mesh from the core
as well as particle fluxes. Returns a dictionary of interpreted results.
"""
function read_b2_boundary_parameters(filename::String)::Dict{String, Any}
    ret_dict = Dict{String, Any}()
    namelist = readnml(filename)
    nbc = namelist[:boundary][:nbc]

    # Sources from core
    ret_dict["power_electrons"] = 0.0  # W
    ret_dict["power_ions"] = 0.0  # W
    ret_dict["number_of_boundaries"] = nbc
    ret_dict["number_of_core_source_boundaries"] = 0

    for bc ∈ 1:nbc
        # Only consider south boundaries. South boundaries are at the interface with
        # the core and at the PFR mesh edges. No power should come from the PFR.
        core_source = false
        if namelist[:boundary][:BCCHAR][bc] == "S"
            # ENE : electron energy condition
            # If BCENE is 8 or 16, ENEPAR's first element will give total power across
            # the boundary in the electron channel.
            # For a double null case, there can be two halves of the core boundary. Just
            # sum all south boundaries and assume someone didn't put power coming in
            # from the PFR.
            # See page 80 of SOLPS user manual 2022 09 13
            bcene = namelist[:boundary][:BCENE][bc]
            bcene_power = [8, 16]
            bcene_pfr_appropriate = [2, 22]
            handled_south_bcene = [bcene_power; bcene_pfr_appropriate]
            if bcene in bcene_power
                enepar = namelist[:boundary][:ENEPAR][bc]
                ret_dict["power_electrons"] += enepar
                if enepar != 0
                    core_source = true
                end
            elseif bcene in bcene_pfr_appropriate
                nothing
            elseif !(bcene in handled_south_bcene)
                throw(
                    ArgumentError(
                        "BCENE type " * repr(bcene) *
                        " is not handled. Acceptable types = " *
                        repr(handled_south_bcene),
                    ),
                )
            end

            # ENI : ion energy condition
            # Similar to ENE, but for ions and there are more options; now 27 should be valid.
            # See page 81 of SOLPS user manual 2022 09 13
            bceni = namelist[:boundary][:BCENI][bc]
            bceni_power = [8, 16, 27]
            bceni_pfr_appropriate = [2, 22]
            handled_south_bceni = [bceni_power; bceni_pfr_appropriate]
            if bceni in bceni_power
                enipar = namelist[:boundary][:ENIPAR][bc]
                ret_dict["power_ions"] += enipar
                if enipar != 0
                    core_source = true
                end
            elseif bceni in bceni_pfr_appropriate
                nothing
            elseif !(bceni in handled_south_bceni)
                throw(
                    ArgumentError(
                        "BCENI type " * repr(bceni) *
                        " is not handled. Acceptable types = " *
                        repr(handled_south_bceni),
                    ),
                )
            end

            #

            if core_source
                ret_dict["number_of_core_source_boundaries"] += 1
            end
        end
    end

    return ret_dict
end
