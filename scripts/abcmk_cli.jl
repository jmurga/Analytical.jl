using Fire,Distributed

"Function to estimate rates"
@main function rates(;ne::Int64=500, samples::Int64=2, workers::Int64=1)
    println("Parsed args:")
    println("opt1=>", ne," opt2=>", samples," workers=>", workers)
    addprocs(2)
    @eval @everywhere using Analytical

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
