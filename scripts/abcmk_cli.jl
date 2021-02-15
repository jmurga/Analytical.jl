using Fire, Distributed

"Function to estimate rates"
@main function rates(;ne::Int64=1000, samples::Int64=500, gamNeg::String="-1000 -200", gL::String="5 10", gH::String="400 1000",dac::String="2,4,5,10,20,50,200,500,700",output::String="/home/jmurga/rates.jld2",workers::Int64=1)

    tmpNeg    = parse.(Int,split(gamNeg," "))
    tmpWeak   = parse.(Int,split(gL," "))
    tmpStrong = parse.(Int,split(gH," "))
    dac       = parse.(Int,split(dac,","))

    addprocs(workers)

    @eval @everywhere using Analytical
    @eval adap = Analytical.parameters(N=$ne,n=$samples,dac=$dac)

    @eval convolutedSamples = Analytical.binomialDict()
    @eval Analytical.binomOp!($adap,$convolutedSamples.bn);
    @time @eval df = Analytical.ratesToStats(param = $adap,convolutedSamples=$convolutedSamples,gH=collect($tmpStrong[1]:$tmpStrong[2]),gL=collect($tmpWeak[1]:$tmpWeak[2]),gamNeg=collect($tmpNeg[1]:$tmpNeg[2]),iterations = 10^3,shape=$adap.al,output=$output);
end

"Summary statistics from rates"
@main function summStat(opt1 = 1, opt2::Int = 2, flag = false)
    println("Parsed args:")
    println("flag=>", flag)
    println("arg=>", x)
    println("opt1=>", opt1)
    println("opt2=>", opt2)
    @eval @everywhere using Analytical
end
