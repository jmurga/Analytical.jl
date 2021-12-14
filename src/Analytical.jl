module Analytical

	using Parameters, SparseArrays, Distributed, CSV, JLD2, DataFrames, ProgressMeter, Quadmath, ParallelUtilities, StatsBase

	# Analytical solutions
	import Roots: find_zero
	import NLsolve: nlsolve
	import SpecialFunctions: polygamma, zeta
	import PoissonRandom: pois_rand
	import Distributions: Binomial, pdf

	# Parse data
	import GZip: open
	import Parsers: parse
	import OrderedCollections: OrderedDict
	import FastaIO: readfasta
	import Random: randstring

	# MK-approaches
	import LsqFit: curve_fit, confidence_interval
	import HypothesisTests: pvalue, FisherExactTest

	include("parameters.jl")
	include("fixations.jl")
	include("polymorphism.jl")
	include("summaryStatistics.jl")
	include("rates.jl")
	include("inferTools.jl")
	include("readFasta.jl")
	include("methods.jl")

end
