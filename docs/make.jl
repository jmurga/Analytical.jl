using Documenter, Analytical


makedocs(
    sitename="Analytical MK approach",
    modules  = [Analytical],
    doctest  = false,
    pages    = [
        "index.md"
    ]
)



deploydocs(
    repo = "github.com/JuliaStats/Distributions.jl.git",
    versions = ["stable" => "v^", "v#.#", "dev" => "master"]
)