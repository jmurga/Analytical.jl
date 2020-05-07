using Documenter, Analytical


makedocs(
    sitename="Analytical MK approach",
    modules  = [Analytical],
    doctest  = false,
    pages    = [
        "index.md"
    ]
)
