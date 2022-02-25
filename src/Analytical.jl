module Analytical


	using Parameters, SparseArrays, Distributed, CSV, JLD2, DataFrames, ProgressMeter, Quadmath, GZip, ParallelUtilities, StatsBase, GSL

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
	import HypothesisTests: pvalue,FisherExactTest
	

	import LinkedLists: LinkedList
	import Printf: @printf
	import QuadGK: quadgk
	import StaticArrays: SVector

	include("parameters.jl")
	include("fixations.jl")
	include("polymorphism.jl")
	include("summary_statistics.jl")
	include("rates.jl")
	include("infer_tools.jl")
	include("read_fasta.jl")
	include("methods.jl")
	include("prefersim.jl")


end
