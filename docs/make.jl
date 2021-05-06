using Documenter,Analytical

makedocs(
	modules   = [Analytical],
	format    = Documenter.HTML(),
	sitename  = "ABC-MK",
	authors   = "Jesús Murga-Moreno, Lawrence Uricchio, David Enard",
	pages     = [
		"Home" => "index.md",
		"Package Overview" => [
			"Analytical estimations" => "analytical.md",
			"Processing data" => "data.md",
			"MK approaches" => "mk.md"
		],
		"Infering the rate and strength of adaptation" =>[
			"Empirical estimation" => "empirical.md",
			"Rates" => "rates.md",
			"Input data" => "input.md",
			"Summary statistics" => "summstat.md",
			"ABC inference" => "abc.md",
			"Maximum-A-Posteriori" => "map.md"
		],	
		"Notebooks" => "notebooks.md",
		"Command line client" => "cli.md",
		"Reference" => "reference.md"
	]
)


deploydocs(
	repo      = "github.com/jmurga/Analytical.jl.git"
)
