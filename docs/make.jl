using Documenter, Analytical

analyticalMakedocs() = makedocs(
	modules = [Analytical],
	sitename="Analytical MK approach",
	authors = "Lawrence Uricchio, Jesús Murga-Moreno",
	format = Documenter.HTML(
		canonical = "https://jmurga.github.io/Analytical.jl/stable/",
	),
	pages = [
		"Home" => "index.md",
		"Package Overview" => [
			"Analytical resolution" => "page.md"
		],
		"References" => [
			"Fixation probabilites" => "reference.md",
		]
	],
)

analyticalMakedocs()

# deploydocs(
#     repo = "github.com/jmurga/Analytical.jl.git"
# )
