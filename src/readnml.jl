"""
This is a bareminimum working example of readnlm and writenlm using package
PyFortran90Namelists as a dependency. This package is registered with Julia Registry
but uses python modules in the backend which we want to avoid.

The codes for readnlm and writenlm have been taken from
https://github.com/ProjectTorreyPines/Fortran90Namelists.jl/blob/d623136e1c91bc0789d675879066b23151d2b74a/src/Namelist.jl
"""

import PyFortran90Namelists: fparse, Tokenizer, lex

fortran_parse(str) = for types ∈ [Int, Float64, Bool, String]
        try
            return fparse(types, str)
        catch
        end
    end

#= ==== =#
#  READ  #
#= ==== =#
"""
    readnml(filename::String; verbose=false)::Dict

Parse fortran namelist in given filename and returns data in nested dictionary structure

NOTE: This parser has the following known limitations (which may be fixed in the future):

  - Cannot handle vector indexes ==> we should use sparsearrays
  - Cannot handle multidimensional arrays
  - Cannot handle complex numbers
  - Cannot handle `!` `;` `#` in strings
  - Cannot handle 1.0+0 exponential notation
  - Will completely neglect comments
  - Will completely neglect text outside of namelist delimiters

These limitations can easily be seen by running regression tests.
Still, even with limited functionalites this should cover most common FORTRAN namelist usage.
"""
function readnml(filename::String; verbose=false)::Dict
    open(filename, "r") do io
        return readnml(io; verbose=verbose)
    end
end

function readnml(io::IO; verbose=false)
    data = Dict()
    return readnml!(io, data; verbose=verbose)
end

function readnml!(io, data; verbose=false)
    tk = Tokenizer()

    h = data

    for line ∈ eachline(io)
        # skip comments or empty lines
        line = split(line, ";")[1]
        line = split(line, "!")[1]
        line = split(line, "#")[1]
        line = strip(line)
        if length(line) == 0
            continue
        end
        line = replace(line, "\$" => "&")
        line = replace(line, r"^&$" => "/")
        line = replace(line, "&end" => "/")

        # remove spaces
        item = [k for k ∈ lex(tk, line) if length(strip(strip(k), ',')) > 0]

        if verbose
            print(strip(line) * " ")
        end

        # open of namelist
        if item[1] == "&"
            if !(item[2] in keys(h))
                h[Symbol(item[2])] = Dict()
                h[Symbol(item[2])][:parent] = h
            end
            h = h[Symbol(item[2])]

            # close of namelist
        elseif item[1] == "/"
            child = h
            h = child[:parent]
            delete!(child, :parent)

            # parsing of elements
        elseif (h !== data)
            if (item[2] == "=")
                # simple values
                if length(item) == 3
                    value = fortran_parse(item[3])

                    # arrays (handles repetitions)
                else
                    tmp = item[3:end]
                    value = Any[]
                    for k ∈ 1:length(tmp)
                        if (k - 1 > 1) && (tmp[k-1] == "*")
                        elseif tmp[k] == "*"
                        elseif (k + 1 < length(tmp)) && (tmp[k+1] == "*")
                            for reps ∈ 1:Int(fortran_parse(tmp[k]))
                                push!(value, fortran_parse(tmp[k+2]))
                            end
                        else
                            push!(value, fortran_parse(tmp[k]))
                        end
                    end
                    value = collect(promote(value...))
                end
                h[Symbol(item[1])] = value
                if verbose
                    print("[$(typeof(value))] -> $(value)")
                end
            else
                if verbose
                    print("[SKIP index]")
                end
            end
        else
            if verbose
                print("[SKIP outside]")
            end
        end
        if verbose
            println()
        end
    end

    return data
end

#= ===== =#
#  WRITE  #
#= ===== =#

"""
    writenml(filename::String, data::Dict; verbose=false)::String

Write nested dictionary structure as fortran namelist to a given filename

NOTE: For a list of known limitations look at the help of readnml()
"""
function writenml(filename::String, data::Dict; verbose=false)::String
    open(filename, "w") do io
        return writenml(io, data; verbose=verbose)
    end
end

function writenml(io::IO, data::Dict; verbose=false)
    txt = []

    for nml ∈ keys(data)
        push!(txt, "&$(nml)")
        for (item, value) ∈ data[nml]
            if typeof(value) <: Vector
                frtn_string =
                    join(map(x -> JuliaToFortran.to_fortran(x).data, value), " ")
            else
                frtn_string = JuliaToFortran.to_fortran(data[nml][item]).data
            end
            push!(txt, "$(item) = $(frtn_string)")
        end
        push!(txt, "/")
    end
    txt = join(txt, "\n")

    if verbose
        println(txt)
    end

    write(io, txt)

    return txt
end