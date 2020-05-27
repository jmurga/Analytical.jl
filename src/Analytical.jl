module Analytical

	include("parameters.jl")
	include("fixations.jl")
	include("polymorphism.jl")
	include("summaryStatistics.jl")
	include("inferTools.jl")
	include("features.jl")

	using Parameters, PyCall, SpecialFunctions, Distributions, Roots, ArbNumerics, StatsBase

	import CSV: read
	import CSV: write
	import DataFrames: DataFrame
	import GZip: open
	import Parsers: parse

	export adap

end
