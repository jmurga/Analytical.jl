using RCall

"""
	Estimating and plotting MAP using locfit and ggplot2 in R. It assume your folder contains the posterior estimated through ABCreg
"""
function plotMap(;analysis::String,output::String)

	try
		@eval using RCall
		@eval R"""library(ggplot2);library(abc)"""

		out = filter(x -> occursin("post",x), readdir(analysis,join=true))
		out = filter(x -> !occursin(".1.",x),out)

		open(x) = Array(CSV.read(GZip.open(x),DataFrame,header=false))
		posteriors = open.(out)

		maxp = DataFrame(Array{Float64,2}(undef,size(posteriors,1),5),[:aw,:as,:a,:gamNeg,:shape])
		R"""getmap <- function(df){
				temp = as.data.frame(df)
			    d <-locfit(~temp[,1],temp);
			    map<-temp[,1][which.max(predict(d,newdata=temp))]
			}"""
		getmap(x) = rcopy(R"""matrix(apply($x,2,getmap),nrow=1)""")
		tmp = getmap.(posteriors)
		maxp = DataFrame(vcat(tmp...),[:aw,:as,:a,:gamNeg,:shape])
	
		al  = maxp[:,1:3]
		gam  = maxp[:,4:end]
		p = rcopy(R"""al = $al
			dal = melt(al)
			pal = ggplot(dal) + geom_density(aes(x=value,fill=variable),alpha=0.5) + scale_fill_manual(values=c('#30504f', '#e2bd9a', '#ab2710'))
			ggsave(pal,filename=paste0($output))
			""")
		return(maxp)
	catch
		println("Please install R, ggplot2 and abc in your system before execute this function")
	end
end
