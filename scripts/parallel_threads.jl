using Parameters, SparseArrays, Distributed, CSV, JLD2, DataFrames, ProgressMeter, Quadmath, GZip, ParallelUtilities, StatsBase, RCall

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

include("parameters.jl")
include("fixations.jl")
include("polymorphism.jl")
include("summaryStatistics.jl")
include("rates.jl")
include("inferTools.jl")
include("readFasta.jl")
include("methods.jl")


function solve(param::parameters,convoluted_samples::binomial_dict,alTot::Float64,alLow::Float64,gH::Int64,gL::Int64,gamNeg::Int64,afac::Float64,θ::Float64,ρ::Float64)
        # Creating model to solve
    # Γ distribution
    param.al = afac; param.be = abs(afac/gamNeg); param.gamNeg = gamNeg
    # α, αW
    param.alLow = alLow; param.alTot = alTot;
    # Positive selection coefficients
    param.gH    = gH;param.gL = gL
    # Mutation rate and recomb
    param.thetaMidNeutral = θ; param.rho = ρ
    #=0.001.thetaMidNeutral = θ; param.θᵣ .= θᵣ; param.rho = ρ=#
    # Solving θ on non-coding region and probabilites to get α value without BGS
    param.B = 0.999
    setThetaF!(param)
    setPpos!(param)    

    return([param.pposL,param.pposH,param.thetaF])

end

function r_old(param::parameters,convoluted_samples::binomial_dict,alTot::Float64,alLow::Float64,gH::Int64,gL::Int64,gamNeg::Int64,afac::Float64,θ::Float64,ρ::Float64)

    # Creating model to solve
    # Γ distribution
    param.al    = afac; param.be = abs(afac/gamNeg); param.gamNeg = gamNeg
    # α, αW
    param.alLow = alLow; param.alTot = alTot;
    # Positive selection coefficients
    param.gH    = gH;param.gL = gL
    # Mutation rate and recomb
    param.thetaMidNeutral = θ; param.rho = ρ
    #=0.001.thetaMidNeutral = θ; param.θᵣ .= θᵣ; param.rho = ρ=#
    # Solving θ on non-coding region and probabilites to get α value without BGS
    param.B = 0.999
    setThetaF!(param)
    setPpos!(param)

    # Allocate array to solve the model for all B values
    r = zeros(size(param.B_bins,1),(size(param.dac,1) * 2) + 14)
    for j in eachindex(param.B_bins)
        # Set B value
        param.B = param.B_bins[j]

        cnvBinom = convoluted_samples.bn[param.B]
        # Solve θ non-coding for the B value.
        setThetaF!(param)
        # Solve model for the B value
        tmp = try
            # Fixation
            fN       = param.B*fixNeut(param)
            fNeg     = param.B*fixNegB(param,0.5*param.pposH+0.5*param.pposL)
            fPosL    = fixPosSim(param,param.gL,0.5*param.pposL)
            fPosH    = fixPosSim(param,param.gH,0.5*param.pposH)

            ds       = fN
            dn       = fNeg + fPosL + fPosH

            # Polymorphism
            neut = DiscSFSNeutDown(param,cnvBinom)
            selH = if isinf(exp(param.gH * 2))
                DiscSFSSelPosDownArb(param,param.gH,param.pposH,cnvBinom)
            else
                DiscSFSSelPosDown(param,param.gH,param.pposH,cnvBinom)
            end
            selL = DiscSFSSelPosDown(param,param.gL,param.pposL,cnvBinom)
            selN = DiscSFSSelNegDown(param,param.pposH+param.pposL,cnvBinom)
            
            # Cumulative rates
            tmp = cumulative_sfs(hcat(neut,selH,selL,selN),false)
            splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));
            neut, selH, selL, selN = splitColumns(tmp)
            sel = (selH+selL)+selN

            ##########
            # Output #
            ##########
           tmp = vcat(param.B,param.alLow,param.alTot,param.gamNeg,param.gL,param.gH,param.al,param.thetaMidNeutral,neut[param.dac],sel[param.dac],ds,dn,fPosL,fPosH,neut[1],sel[1])'
        catch e
            tmp = zeros(32)
        end
        @inbounds r[j,:] = tmp
        #=@inbounds r = tmp=#
    end

    return r
end

function r(param::parameters,convoluted_samples::binomial_dict,alTot::Float64,alLow::Float64,gH::Int64,gL::Int64,gamNeg::Int64,afac::Float64,θ::Float64,ρ::Float64,e::Vector{Float64})

    # Creating model to solve
    # Γ distribution
    param.al    = afac; param.be = abs(afac/gamNeg); param.gamNeg = gamNeg
    # α, αW
    param.alLow = alLow; param.alTot = alTot;
    # Positive selection coefficients
    param.gH    = gH;param.gL = gL
    # Mutation rate and recomb
    param.thetaF = e[3]; param.thetaMidNeutral = θ; param.rho = ρ
    # Fixation probabilites
    param.pposL = e[1]; param.pposH = e[2]

    # Solving θ on non-coding region and probabilites to get α value without BGS
    param.B = e[4]
    cnvBinom = convoluted_samples.bn[e[4]]
    # Solve θ non-coding for the B value.
    setThetaF!(param)
    # Solve model for the B value

    # Fixation
    fN       = param.B*fixNeut(param)
    fNeg     = param.B*fixNegB(param,0.5*param.pposH+0.5*param.pposL)
    fPosL    = fixPosSim(param,param.gL,0.5*param.pposL)
    fPosH    = fixPosSim(param,param.gH,0.5*param.pposH)

    ds       = fN
    dn       = fNeg + fPosL + fPosH

    # Polymorphism
    neut = DiscSFSNeutDown(param,cnvBinom)
    selH = if isinf(exp(param.gH * 2))
        DiscSFSSelPosDownArb(param,param.gH,param.pposH,cnvBinom)
    else
        DiscSFSSelPosDown(param,param.gH,param.pposH,cnvBinom)
    end
    selL = DiscSFSSelPosDown(param,param.gL,param.pposL,cnvBinom)
    selN = DiscSFSSelNegDown(param,param.pposH+param.pposL,cnvBinom)

    # Cumulative rates
    tmp = cumulative_sfs(hcat(neut,selH,selL,selN),false)
    splitColumns(matrix::Array{Float64,2}) = (view(matrix, :, i) for i in 1:size(matrix, 2));
    neut, selH, selL, selN = splitColumns(tmp)
    sel = (selH+selL)+selN

    ##########
    # Output #
    ##########
    r = vcat(param.B,param.alLow,param.alTot,param.gamNeg,param.gL,param.gH,param.al,param.thetaMidNeutral,neut[param.dac],sel[param.dac],ds,dn,fPosL,fPosH,neut[1],sel[1])'
    
    return r
end

param = parameters()
convoluted_samples = binomial_dict()
binomOp!(param,convoluted_samples.bn)

#gH=200:1000;gL=1:10;gamNeg=-2000:-200;iterations = 10^2;shape=adap.al;theta=0.001;rho=0.001

iterations = 10^3
afac       = fill(param.al,iterations)

# Random α values
nTot       = fill(0.4,iterations)
# Random weak selection coefficients
nLow       = fill(0.2,iterations)
# Creating N models to iter in threads. Set N models (paramerters) and sampling probabilites (binomial_dict)
nParam     = [param for i in 1:iterations];
nBinom     = [convoluted_samples for i in 1:iterations];
# Random strong selection coefficients
ngh        = fill(200,iterations)
ngl        = fill(5,iterations)
# Random negative selection coefficients
ngamNeg    = fill(-457,iterations)
θ          = fill(0.001,iterations)
ρ          = fill(0.001,iterations)

function rates_threads(;param::parameters,
                convoluted_samples::binomial_dict,
                gH::S,
                gL::S,
                gamNeg::S,
                theta::Union{Float64,Nothing}=0.001,
                rho::Union{Float64,Nothing}=0.001,
                shape::Float64=0.184,
                iterations::Int64,
                output::String) where S <: Union{Array{Int64,1},UnitRange{Int64},Nothing}
    

    # Iterations = models to solve
    # Factor to modify input Γ(shape) parameter. Flexible Γ distribution over negative alleles
    fac     = rand(-2:0.05:2,iterations)
    afac    = @. param.al*(2^fac)
    
    # Deleting shape > 1. Negative alpha_x values
    idx = findall(afac .> 1)
    if !isempty(idx)
        afac[idx] = rand(afac[afac .< 1],size(idx,1))
    end

    # Random α values
    nTot    = rand(0.1:0.01:0.9,iterations)
    
    # Defining αW. It is possible to solve non-accounting for weak fixations
    if isnothing(gL)
        # Setting αW to 0 for all estimations
        nLow    = fill(0.0,iterations)
        # Random strong selection coefficients
        ngl     = rand(repeat([1],iterations),iterations);
    else
        # Setting αW as proportion of α
        lfac    = rand(0.0:0.05:0.9,iterations)
        nLow    = @. nTot * lfac
        # Random weak selection coefficients
        ngl     = rand(repeat(gL,iterations),iterations);
    end

    # Creating N models to iter in threads. Set N models (paramerters) and sampling probabilites (binomial_dict)
    nParam  = [param for i in 1:iterations];
    nBinom  = [convoluted_samples for i in 1:iterations];
    
    # Random strong selection coefficients
    ngh     = rand(repeat(gH,iterations),iterations);
    # Random negative selection coefficients
    ngamNeg = rand(repeat(gamNeg,iterations),iterations);

    # Random θ on coding regions
    #=obsNeut        = sfs.p0./sum(sfs.p0) 
    n              = 1 ./collect(1:(param.nn-1))
    expNeut        = n ./ sum(n)=#

    if !isnothing(theta)
        θ = fill(theta,iterations)
        #=θᵣ= fill(theta .* obsNeut ./ expNeut,iterations)=#
    else
        θ = rand(0.0005:0.0005:0.01,iterations)
        #=θᵣ= map(x -> x .* obsNeut ./ expNeut,θ)=#
    end

    # Random ρ on coding regions
    if !isnothing(rho)
        ρ = fill(rho,iterations)
    else
        ρ = rand(0.0005:0.0005:0.05,iterations)
    end

    ######Parallel; 10^5 -> 62.76"
    x = zeros(1,3,iterations);
    @inbounds @sync for i in 1:iterations
        Base.Threads.@spawn begin
            @inbounds x[:,:,i] = solve(param,convoluted_samples, nTot[i], nLow[i], ngh[i], ngl[i], ngamNeg[i], afac[i], θ[i], ρ[i])
        end
    end    
    
    #y = ThreadsX.map(solve, nParam,nBinom, nTot, nLow, ngh, ngl, ngamNeg, afac, θ, ρ)

    x = repeat(x, outer = [size(param.B_bins,1), 1, 1]);
    x = vcat([ @view x[:,:,i] for i=1:iterations]...);
    b = repeat(param.B_bins,iterations);
    e = hcat(x,b);

    out = zeros(iterations,32,size(param.B_bins,1));

    @showprogress for (j,val) in enumerate(reverse(param.B_bins))
        tmp = e[e[:,4] .== val,:]
        for i in 1:iterations
            Base.Threads.@spawn begin
               out[i,:,j] = r(nParam[i],nBinom[i], nTot[i], nLow[i], ngh[i], ngl[i], ngamNeg[i], afac[i], θ[i], ρ[i],tmp[i,:]);
            end
        end
    end

    df = vcat([@view out[:,:,i] for i=1:size(out,3)]...);

    # Saving models and rates
    models = DataFrame(df[:,1:8],[:B,:alLow,:alTot,:gamNeg,:gL,:gH,:al,:ρ]);
    neut   = df[:,9:(8+size(param.dac,1))];
    sel    = df[:,(9+size(param.dac,1)):(8+size(param.dac,1)*2)];
    dsdn   = Array(df[:,(end-5):end-2]);

    sum_pol = df[:,end-1:end];

    # Saving multiple summary statistics
    n = OrderedDict{Int,Array}();
    s = OrderedDict{Int,Array}();
    for i in eachindex(param.dac)
        n[param.dac[i]] = neut[:,i]
        s[param.dac[i]] = sel[:,i]
    end

    #=return models=#
    # Writting HDF5 file
    JLD2.jldopen(output, "a+") do file
        file[string(param.N)* "/" * string(param.n) * "/models"] = models;
        file[string(param.N)* "/" * string(param.n) * "/neut"]   = n;
        file[string(param.N)* "/" * string(param.n) * "/sel"]    = s;
        file[string(param.N)* "/" * string(param.n) * "/dsdn"]   = dsdn;
        file[string(param.N)* "/" * string(param.n) * "/pol"]   = sum_pol;
        file[string(param.N)* "/" * string(param.n) * "/dac"]    = param.dac;
    end;
end

@time df = rates_threads(param = adap,convoluted_samples=convoluted_samples,gH=200:2000,gL=1:10,gamNeg=-2000:-200,iterations = 10^3,shape=adap.al,output="/home/jmurga/rates_threads.jld2");rm("/home/jmurga/rates_threads.jld2")

@time df = Analytical.rates(param = adap,convoluted_samples=convoluted_samples,gH=200:2000,gL=1:10,gamNeg=-2000:-200,iterations = 10^2,shape=adap.al,output="/home/jmurga/rates_distributed.jld2");rm("/home/jmurga/rates_distributed.jld2")

######Sequential

real = fill(Matrix{Float64}(undef, 0,0),iterations);
@time for i in eachindex(nTot)
    real[i] = r(param,convoluted_samples, nTot[i], nLow[i], ngh[i], ngl[i], ngamNeg[i], afac[i], θ[i], ρ[i])
end


### ThreadsPools
#=@time @sync begin @qbthreads for i in 1:iterations
        x[:,:,i] = solve(nParam[i],nBinom[i], nTot[i], nLow[i], ngh[i], ngl[i], ngamNeg[i], afac[i], θ[i], ρ[i])
    end
end=#
#=@showprogress for (j,val) in enumerate(param.B_bins)

    tmp = e[e[:,4] .== val,:]
    @sync begin @qbthreads for i in 1:iterations
           out[i,:,j] = r(nParam[i],nBinom[i], nTot[i], nLow[i], ngh[i], ngl[i], ngamNeg[i], afac[i], θ[i], ρ[i],tmp[i,:])
        end
    end
end
=#

