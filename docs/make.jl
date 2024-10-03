using Documenter
using SOLPS2imas

makedocs(;
    modules=[SOLPS2imas],
    format=Documenter.HTML(),
    sitename="SOLPS2imas",
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/ProjectTorreyPines/SOLPS2imas.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="master",
    versions=["stable" => "v^", "v#.#"],
)
