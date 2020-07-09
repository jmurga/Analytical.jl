using Analytical, CSV, DataFrames 

function analyticalApproach(param)

    B = param.B
	Analytical.set_theta_f(param)
	theta_f = param.theta_f
	param.B = 0.999
	Analytical.set_theta_f(param)
	Analytical.setPpos(param)
	param.theta_f = theta_f
	param.B = B
	x,y = Analytical.analyticalAlpha(param=param)

	return round.(hcat(theta_f/(4*param.N),param.pposL,param.pposH,param.alLow,param.alTot,round(y[end],digits=5),param.B),digits=12)
end


alpha = [ 0.2 , 0.3, 0.4]
bgs   = [0.2,0.4,0.8,0.999]

function simTable(alphas,bgsValues,pSize)

    out = zeros(size(alphas,1)*size(bgsValues,1)*3,7)
    it = 1

    for a in 1:size(alphas,1)
        for b in 1:size(bgsValues,1)

            adap1 = Analytical.parameters(N=pSize,n=661,gam_neg=-457,Lf=2*10^5, B=bgsValues[b],gL=10,gH=500,alTot=alphas[a],alLow=alphas[a]*0.25)

            adap2 = Analytical.parameters(N=pSize,n=661,gam_neg=-457,Lf=2*10^5, B=bgsValues[b],gL=10,gH=500,alTot=alphas[a],alLow=alphas[a]*0.5)

            adap3 = Analytical.parameters(N=pSize,n=661,gam_neg=-457,Lf=2*10^5, B=bgsValues[b],gL=10,gH=500,alTot=alphas[a],alLow=alphas[a]*0.75)

            r1 = analyticalApproach(adap1)
            r2 = analyticalApproach(adap2)
            r3 = analyticalApproach(adap3)
            
            out[it,:] = r1
            out[it+1,:] = r2
            out[it+2,:] = r3
            it = it + 3
      end
    end
    out = sort(DataFrame(out),5)
    out[:,4] = round.(out[:,4],digits=2)

    rename!(out, [Symbol("bgsThetaF"),Symbol("pposL"),Symbol("pposH") ,Symbol("alphaW"), Symbol("alpha"),Symbol("estimation"),Symbol("B")])
    
    return out
end    

simulations = simTable(alpha,bgs,parse(Int,ARGS[1]))
println(parse(Int,ARGS[1]))
CSV.write("/home/jmurga/mkt/202004/rawData/simulations/simTable.tsv", simulations, delim='\t')
