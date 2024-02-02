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
            key = strip(split(line)[1], ['\''])
            value = split(line, '\'', keepempty=false)[2:end]
            try
                value = [
                    '.' in v ? parse(Float64, v) : parse(Int, v)
                    for v in value if length(strip(v, ' ')) > 0
                ]
            catch
                # Adapted from the method for parsing b2mn or b2ag in omfit_solps.py
                # I don't know when or why this is needed, but that parser has been
                # tested aggressively and works well, so I'm copying it.
                value = split(line, '\'')[4]
            end
            value = length(value) == 1 ? value[1] : value
            contents[key] = value
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
