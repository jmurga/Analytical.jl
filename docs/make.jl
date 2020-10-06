using Documenter, AnalyticaluSfsFromFasta


makedocs(
    modules   = [Analytical],
    format    = Documenter.HTML(),
    sitename  = "ABC-MK",
    authors   = "Jesús Murga-Moreno, Lawrence Uricchio, David Enard",
    pages     = [
        "Home" => "index.md",
        "Package Overview" => [
            "Analytical estimations" => "analytical.md",
            "ABC inference" => "abc.md"
            "Raw data" => "data.md"
        ],
        "Reference" => "reference.md"
    ]
)


deploydocs(
    repo      = "github.com/jmurga/Analytical.jl.git"
)
