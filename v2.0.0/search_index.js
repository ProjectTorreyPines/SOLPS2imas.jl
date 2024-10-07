var documenterSearchIndex = {"docs":
[{"location":"#SOLPS2imas.jl","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"","category":"section"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"Pages = [\"index.md\"]\nDepth = 5","category":"page"},{"location":"#Installation","page":"SOLPS2imas.jl","title":"Installation","text":"","category":"section"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"For installation:","category":"page"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"using Pkg\nPkg.add(\"SOLPS2imas\")","category":"page"},{"location":"#solps2imas","page":"SOLPS2imas.jl","title":"solps2imas","text":"","category":"section"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"The main function of this module is solps2imas which reads SOLPS input and output files and creates an IMAS IDS instance.","category":"page"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"solps2imas","category":"page"},{"location":"#SOLPS2imas.solps2imas","page":"SOLPS2imas.jl","title":"SOLPS2imas.solps2imas","text":"solps2imas(\n    b2gmtry::String,\n    b2output::String=\"\";\n    gsdesc::String=default_gsdesc,\n    b2mn::String=\"\",\n    fort::Tuple{String, String, String}=(\"\", \"\", \"\"),\n    fort_tol::Float64=1e-6,\n    b2_parameters::Tuple{String, String, String, String,}=(\"\", \"\", \"\", \"\"),\n    load_bb::Bool=false,\n)::IMASdd.dd\n\nMain function of the module. Takes in a geometry file and an (optional) output file (either b2time or b2fstate) and a grid_ggd description in the form of a Dict or filename to equivalent YAML file. Additionally, EIRENE fort files can be provided as tuple of 3 filenames consisting fort.33, fort.34, and fort.35 files. The grids in these files are matched with SOLPS grid with a tolerance of fort_tol (defaults to 1e-6). Further settings can be loaded from b2.*.parameters files and equilibrium files.\n\n\n\n\n\n","category":"function"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"The main solps2imas feature can be called to load everything. However, some special cases may demand only subsets of the data. In these cases, it is possible to call subcomponents of the SOLPS to IMAS data transfer feature by using these functions:","category":"page"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"load_summary_data!","category":"page"},{"location":"#SOLPS2imas.load_summary_data!","page":"SOLPS2imas.jl","title":"SOLPS2imas.load_summary_data!","text":"load_summary_data!(\n    ids::IMASdd.dd,\n    b2_parameters::Tuple{String, String, String, String}=(\"\", \"\", \"\", \"\");\n)\n\nLoads high level shot summary data into the summary IDS after reading and interpreting SOLPS input files, such as b2.*.parameters.\n\n\n\n\n\n","category":"function"},{"location":"#Parsing-tools","page":"SOLPS2imas.jl","title":"Parsing tools","text":"","category":"section"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"This module uses some parsing functions which can be used standalone as well","category":"page"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"read_b2_output\nread_b2mn_output\nread_b2time_output","category":"page"},{"location":"#SOLPS2imas.read_b2_output","page":"SOLPS2imas.jl","title":"SOLPS2imas.read_b2_output","text":"read_b2_output(filename::String)::Dict{String, Dict{String, Any}}\n\nRead final state b2 output file (b2fstate or b2time.nc) or b2fgmtry file and return a dictionary with structure: Dict(\"dim\" => Dict{String, Any}, \"data\" => Dict{String, Any}) where \"dim\" contains the dimensions of the data and \"data\" contains the data itself, with keys corresponding to the field names.\n\nSupported SOLPS files as input via filename:\n\nb2fstate\nb2fstati\nb2time.nc\nb2fgmtry\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.read_b2mn_output","page":"SOLPS2imas.jl","title":"SOLPS2imas.read_b2mn_output","text":"read_b2mn_output(filename::String)::Dict{String, Any}\n\nRead b2mn output file and store the quantities in a dictionary.\n\nSupported SOLPS files as input via filename:\n\nb2mn.dat\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.read_b2time_output","page":"SOLPS2imas.jl","title":"SOLPS2imas.read_b2time_output","text":"read_b2time_output(filename::String)::Dict{String, Dict{String, Any}}\n\nRead time dependent b2 output file and return a dictionary with structure: Dict(\"dim\" => Dict{String, Any}, \"data\" => Dict{String, Any}) where \"dim\" contains the dimensions of the data and \"data\" contains the data itself, with keys corresponding to the field names.\n\nSupported SOLPS files as input via filename:\n\nb2time.nc\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS-Mesh-Tools","page":"SOLPS2imas.jl","title":"SOLPS Mesh Tools","text":"","category":"section"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"SOLPS2imas.xytoc\nSOLPS2imas.ctoxy\nSOLPS2imas.data_xytoc\nSOLPS2imas.search_point\nSOLPS2imas.search_edges\nSOLPS2imas.distance_between_nodes\nSOLPS2imas.neighbour_inds\nSOLPS2imas.get_neighbour_inds\nSOLPS2imas.attach_neightbours!","category":"page"},{"location":"#SOLPS2imas.xytoc","page":"SOLPS2imas.jl","title":"SOLPS2imas.xytoc","text":"xytoc(ix::Int, iy::Int; nx::Int)::Int\n\nConverts SOLPS indices for crx, cry (ix, iy) that go from 1:nx, 1:ny into the linear index ic used in IMAS for corresponding cells.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.ctoxy","page":"SOLPS2imas.jl","title":"SOLPS2imas.ctoxy","text":"ctoxy(ic::Int; nx::Int)::Tuple{Int, Int}\n\nInverse of xytoc.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.data_xytoc","page":"SOLPS2imas.jl","title":"SOLPS2imas.data_xytoc","text":"data_xytoc(data::Matrix{Float64}; nx::Int)::Vector{Float64}\n\nFlattens 2d data given on cell indices ix, iy into 1d data on linear index ic. ic is calculated using xytoc function. Data is assumed to have dimensions (ny, nx) where ny is not required in this conversion.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.search_point","page":"SOLPS2imas.jl","title":"SOLPS2imas.search_point","text":"search_point(\n    nodes::IMASdd.IDSvector{IMASdd.edge_profiles__grid_ggd___space___objects_per_dimension___object{T}},\n    r::Real,\n    z::Real;\n    tol::Float64=0.0,\n)::Int where {T}\n\nSearch if a point (r, z) is present in the nodes array. Here nodes is generally available in ids.edge_profiles.grid_ggd[:].space[:].objects_per_dimension[1].object\n\nIf the point is not found, the function returns 0.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.search_edges","page":"SOLPS2imas.jl","title":"SOLPS2imas.search_edges","text":"search_edges(\n    edges::IMASdd.IDSvector{IMASdd.edge_profiles__grid_ggd___space___objects_per_dimension___object{T}},\n    edge_nodes::Array{Int, 1}\n)::Int where {T}\n\nSearch if an edge with nodes as edge_nodes already exists in the edges array. edges is generally available in ids.edge_profiles.grid_ggd[:].space[:].objects_per_dimension[2].object\n\nIf the edge is not found, the function returns 0.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.distance_between_nodes","page":"SOLPS2imas.jl","title":"SOLPS2imas.distance_between_nodes","text":"distance_between_nodes(nodes::edges_nodes_type, node_inds::Array{Int, 1})\n\nReturn distance between two nodes with indices node_inds in nodes array. nodes is generally available in ids.edge_profiles.grid_ggd[:].space[:].objects_per_dimension[1].object.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.neighbour_inds","page":"SOLPS2imas.jl","title":"SOLPS2imas.neighbour_inds","text":"neighbour_inds(\n    ic::Int;\n    nx::Int,\n    ny::Int,\n    leftcut::Int,\n    rightcut::Int,\n    topcut::Int,\n    bottomcut::Int,\n)::Vector{Int}\n\n(deprecated function)\n\nReturns indices of neighbours of cell with linear index ic. This function uses the SOLPS grid generation algorithm to determine the neighbours. However, SOLPS geometry file actually provides the neighbor indices directly. Thus, this function is not used in the code anywhere but is kept here for reference.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.get_neighbour_inds","page":"SOLPS2imas.jl","title":"SOLPS2imas.get_neighbour_inds","text":"get_neighbour_inds(\n    ic::Int,\n    gmtry::Dict{String, Dict{String, Any}},\n    it::Int,\n)::Vector{Int}\n\nReturns indices of neighbours of cell with linear index ic. This function uses the parsed SOLPS geometry file to determine the neighbours by using matrices named as leftix, rightix, topix, bottomix, leftiy, rightiy, topiy, and bottomiy.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.attach_neightbours!","page":"SOLPS2imas.jl","title":"SOLPS2imas.attach_neightbours!","text":"attach_neightbours!(\n    cells::IMASdd.IDSvector{IMASdd.edge_profiles__grid_ggd___space___objects_per_dimension___object{T}},\n    edges::IMASdd.IDSvector{IMASdd.edge_profiles__grid_ggd___space___objects_per_dimension___object{T}},\n    gmtry::Dict{String, Dict{String, Any}},\n    it::Int,\n) where {T}\n\nThis function attaches neighbours to each boundary of each cell and each boundary of each edge using the parsed SOLPS geometry file.\n\n\n\n\n\n","category":"function"},{"location":"#Subset-Identification-Tools","page":"SOLPS2imas.jl","title":"Subset Identification Tools","text":"","category":"section"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"These functions provide a way to identify if a linear cell index in SOLPS notation belongs to a grid subset.","category":"page"},{"location":"","page":"SOLPS2imas.jl","title":"SOLPS2imas.jl","text":"SOLPS2imas.in_core\nSOLPS2imas.in_sol\nSOLPS2imas.in_idr\nSOLPS2imas.in_odr\nSOLPS2imas.is_x_aligned\nSOLPS2imas.is_y_aligned\nSOLPS2imas.is_core_cut\nSOLPS2imas.is_pfr_cut\nSOLPS2imas.is_outer_throat\nSOLPS2imas.is_inner_throat\nSOLPS2imas.is_outer_midplane\nSOLPS2imas.is_inner_midplane\nSOLPS2imas.is_outer_target\nSOLPS2imas.is_inner_target\nSOLPS2imas.is_core_boundary\nSOLPS2imas.is_separatrix\nSOLPS2imas.get_xpoint_nodes","category":"page"},{"location":"#SOLPS2imas.in_core","page":"SOLPS2imas.jl","title":"SOLPS2imas.in_core","text":"in_core(; ix, iy, topcut, bottomcut, leftcut, rightcut)::Bool\n\nReturns true if cell indexed ix, iy lie inside the core.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.in_sol","page":"SOLPS2imas.jl","title":"SOLPS2imas.in_sol","text":"in_sol(; iy, topcut, kwargs...)::Bool\n\nReturns true if cell indexed ix, iy lie inside the SOL.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.in_idr","page":"SOLPS2imas.jl","title":"SOLPS2imas.in_idr","text":"in_idr(; ix, iy, topcut, bottomcut, leftcut, kwargs...)::Bool\n\nReturns true if cell indexed ix, iy lie inside the inner divertor region.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.in_odr","page":"SOLPS2imas.jl","title":"SOLPS2imas.in_odr","text":"in_odr(; ix, iy, topcut, bottomcut, rightcut, kwargs...)::Bool\n\nReturns true if cell indexed ix, iy lie inside the outer divertor region.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_x_aligned","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_x_aligned","text":"is_x_aligned(;boundary_ind)::Bool\n\nx_aligned edges will have odd boundary_ind based on chosen order of numbering them.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_y_aligned","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_y_aligned","text":"is_y_aligned(; boundary_ind)::Bool\n\ny_aligned edges will have even boundary_ind based on chosen order of numbering them.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_core_cut","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_core_cut","text":"is_core_cut(;\n    ix,\n    iy,\n    cells,\n    nx,\n    boundary_ind,\n    topcut,\n    bottomcut,\n    leftcut,\n    rightcut,\n)::Bool\n\nReturns true if boundary_ind of a cell at ix, iy is on core_cut (Y-aliged edge).\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_pfr_cut","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_pfr_cut","text":"is_pfr_cut(;\n    ix,\n    iy,\n    cells,\n    nx,\n    boundary_ind,\n    topcut,\n    bottomcut,\n    leftcut,\n    rightcut,\n)::Bool\n\nReturns true if boundary_ind of a cell at ix, iy is on core_cut (y-aliged edge).\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_outer_throat","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_outer_throat","text":"is_outer_throat(; ix, iy, boundary_ind, topcut, rightcut, kwargs...)::Bool\n\nReturns true if boundary_ind of a cell at ix, iy is on outer throat.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_inner_throat","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_inner_throat","text":"is_inner_throat(; ix, iy, boundary_ind, topcut, leftcut, kwargs...)::Bool\n\nReturns true if boundary_ind of a cell at ix, iy is on outer throat.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_outer_midplane","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_outer_midplane","text":"is_outer_midplane(; ix, jxa, boundary_ind)\n\nReturns true if boundary_ind of a cell at ix, iy is on outer midplane.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_inner_midplane","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_inner_midplane","text":"is_inner_midplane(; ix, jxa, boundary_ind)\n\nReturns true if boundary_ind of a cell at ix, iy is on outer midplane.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_outer_target","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_outer_target","text":"is_outer_target(; ix, nx, boundary_ind)\n\nReturns true if boundary_ind of a cell at ix, iy is on outer target.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_inner_target","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_inner_target","text":"is_inner_target(; ix, boundary_ind, kwargs...)::Bool\n\nReturns true if boundary_ind of a cell at ix, iy is on inner target.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_core_boundary","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_core_boundary","text":"is_core_boundary(;\nix,\niy,\nboundary_ind,\nbottomcut,\nleftcut,\nrightcut,\nkwargs...,\n\n)\n\nReturns true if boundary_ind of a cell at ix, iy is on core boundary (central blank spot boundary).\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.is_separatrix","page":"SOLPS2imas.jl","title":"SOLPS2imas.is_separatrix","text":"is_separatrix(; iy, boundary_ind, topcut, kwargs...)::Bool\n\nReturns true if boundary_ind of a cell at ix, iy is on separatrix.\n\n\n\n\n\n","category":"function"},{"location":"#SOLPS2imas.get_xpoint_nodes","page":"SOLPS2imas.jl","title":"SOLPS2imas.get_xpoint_nodes","text":"get_xpoint_nodes(\n    gmtry::Dict{String, Dict{String, Any}},\n)::Vector{Vector{Vector{Float64}}}\n\nLimited to finding first x-point for now. Returns x-point (r, z) for each time index of grid_ggd for the first x-point only. Thus second index correspond to the rank of x-point which is always 1 from output of this function for now.\n\n\n\n\n\n","category":"function"}]
}
