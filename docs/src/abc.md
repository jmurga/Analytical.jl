# *ABC* inference from empirical data

## Prior distributions from analytical estimations

In this example we provide a solution to replicate results at [Uricchio et al. 2019](https://doi.org/10.1038/s41559-019-0890-6). We will simulate $10^6$ summary statistics from random *DFE* to use as prior distribution in *ABCreg*. In this case we will need a set of empirical observed values in order to subset the summary statistics.

We need to set the model accounting for the sampling value. The SFS is expected to be in raw frequencies. If the model is not properly set up, the SFS will not be correctly parsed. In our case, we are going to set up a model with default parameters only to parse the SFS and convolute the observed frequencies with a binomial distribution.

```julia
Analytical.changeParameters(N=1000,n=661,convoluteBinomial=true)
```

Once the model account for the number of samples we can open the files. The function `Analytical.parseSfs` will return polymorphic and divergent counts and SFS accounting for the whole spectrum: `collect(1:adap.nn)/adap.nn`. In addition an output file will be created contained the observed values to input in *ABCreg*.

```julia
path= "/home/jmurga/mktest/data/";suffix="txt";
files = path .* filter(x -> occursin(suffix,x), readdir(path))

empiricalValues = Analytical.parseSfs(data=files,output="testData.tsv",sfsColumns=[3,5],divColumns=[6,7])
```

We make a function to perform $10^6$ simulated values. We solve the analytical approximation taking into account totally random and independent to draw *DFE* and $\alpha_{(x)}$. Each parameters combination are replicated to 5% frequency bins background selection values (saved at `adap.bRange`). 

In Julia you can easily parallelize a loop using ```$ export JULIA_NUM_THREADS=8```. Each iteration will be executed in a thread. In order to check the threads configured, just use in the julia console ```julia> Threads.nthreads()``` before the execution. We compute this example in a Intel i7-7700HQ (8) @ 3.800GHz laptop with 16GB of RAM using 8 threads. Please check [parallelization manual](https://docs.julialang.org/en/v1/manual/parallel-computing/) in order to send the process in a multicore system (or just put two process manually the [*alphaSumStats.jl*](https://github.com/jmurga/Analytical.jl/blob/master/scripts/alphaSumStats.jl) , a script provided to launch from command line).

```julia
function summStats(iter,data)
#@threads
	@showprogress for i in 1:iter
		gam_neg   = -rand(80:400)
		gL        = rand(10:20)
		gH        = rand(200:500)
		alLow     = rand(collect(0.0:0.1:0.2))
		alTot     = rand(collect(0.0:0.05:0.2))

		for j in adap.bRange
			Analytical.changeParameters(
				gam_neg             = gam_neg,
				gL                  = gL,
				gH                  = gH,
				alLow               = alLow,
				alTot               = alTot,
				theta_f             = 1e-3,
				theta_mid_neutral   = 1e-3,
				al                  = 0.184,
				be                  = 0.000402,
				B                   = j,
				bRange              = adap.bRange,
				pposL               = 0.001,
				pposH               = 0.0,
				N                   = 1000,
				n                   = 661,
				Lf                  = 10^6,
				rho                 = 0.001,
				TE                  = 5.0,
				convoluteBinomial   = false
			)

			Analytical.set_theta_f()
			theta_f        = adap.theta_f
			adap.B         = 0.999
			Analytical.set_theta_f()
			Analytical.setPpos()
			adap.theta_f   = theta_f
			adap.B         = j

			x,y,z          = Analytical.alphaByFrequencies(gammaL=adap.gL,gammaH=adap.gH,pposL=adap.pposL,pposH=adap.pposH,observedData=data)
			Analytical.summaryStatistics("/home/jmurga/prior.csv", z)
		end
	end
end

# Launch 10^6 solutions
@time summStats(58824,empiricalValues)
```

## *ABC* inference
Generic ABC methods proceed by three main steps: (1) first sampling parameter values from prior distributions, (2) next simulating a model and calculating informative summary statistics, and lastly (3) comparing the simulated summary statistics to observed data. The parameter values that produce summary statistics that best match the observed data form an approximate posterior distribution. We link Julia with *ABCreg*. It will output one file per line in data. The files contain the posterior distributions. We return the posterior distributions, mean and quantiles.

```julia
posteriors, results  = Analytical.ABCreg(data="/home/jmurga/data.tsv",prior="/home/jmurga/prior.tsv", nparams=27, nsummaries=24, outputPath="/home/jmurga/", outputPrefix="outPaper", tolerance=0.001, regressionMode="T",regPath="/home/jmurga/ABCreg/src/reg")
```

You can easily plot the posterior distributions using Julia or just input the files at your favorite plot software.

```julia
using Plots, Plots.PlotMeasures, StatsPlots, LaTeXStrings

function plotPosterior(data,file,imgSize)

	Plots.gr()
	Plots.theme(:wong2)

	p1 = StatsPlots.density(data],
							legend    = :topright,
							fill      = (0, 0.3),
							xlabel    = L"\alpha",
							label     = [L"\alpha_S" L"\alpha_W" L"\alpha"],
							ylabel    = "Posterior density",
							lw        = 0.5,
							fmt       = :svg,
							size      = imgSize
						)
	Plots.savefig(file)

	return p1
end

p = plotPosterior(posterior[1],"/home/jmurga/posterior1.svg",(800,600))
```

![image](https://raw.githubusercontent.com/jmurga/Analytical.jl/master/docs/src/fig2.svg)
