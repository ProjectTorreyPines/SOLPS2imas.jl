using Documenter
using SOLPS2IMAS

makedocs(;
    modules=[SOLPS2IMAS],
    format=Documenter.HTML(),
    sitename="SOLPS2IMAS",
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/ProjectTorreyPines/SOLPS2IMAS.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="dev",
    versions=["stable" => "v^", "v#.#"],
)
