module Analytical

include("parameters.jl")
include("fixations.jl")
include("polymorphism.jl")
include("summaryStatistics.jl")
include("inferTools.jl")
include("features.jl")

using Parameters, PyCall, SpecialFunctions, Distributions, Roots, CSV, Parsers, StatsBase, Plots, StatsPlots
import Plots.PlotMeasures

export adap

end # module

