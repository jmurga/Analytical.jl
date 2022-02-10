try
	using RCall
	R"""library(locfit);library(ggplot2)"""
catch
	using Pkg
	ENV["R_HOME"]="*"

	Pkg.add("Conda")

	using Conda
	Conda.add("r-base",channel="conda-forge")
	Conda.add(["r-locfit","r-ggplot2","r-data.table","r-r.utils"],channel="conda-forge")

	Pkg.add("RCall")

	using RCall
end

"""
	Estimating and plotting MAP using locfit and ggplot2 in R. It assume your folder contains the posterior estimated through ABCreg
"""
function plot_map(;analysis_folder::String,weak::Bool=true,title::String="Posteriors")

	R"""library(locfit);library(ggplot2);library(data.table)""";

	out = filter(x -> occursin("post",x), readdir(analysis_folder,join=true))
	out          = filter(x -> !occursin(".1.",x),out)

	open(x)      = Array(read(x,DataFrame,header=false))

	# Control outlier inference. 2Nes non negative values
	flt(x)       = x[x[:,4] .> 0,:]
	posteriors   = flt.(open.(out))

	R"""getmap <- function(df){
			temp = as.data.frame(df)
			d <-locfit(~temp[,1],temp);
			map<-temp[,1][which.max(predict(d,newdata=temp))]
		}"""

	getmap(x)    = rcopy(R"""suppressWarnings(matrix(apply($x,2,getmap),nrow=1))""")
	
	if !weak
		posteriors = [posteriors[i][:,3:end] for i in eachindex(posteriors)]
		tmp          = getmap.(posteriors)
		maxp         = DataFrame(vcat(tmp...),[:a,:gamNeg,:shape])
		al           = maxp[:,1:1]
		gam          = maxp[:,2:end]
	else
		tmp          = getmap.(posteriors)
		maxp         = DataFrame(vcat(tmp...),[:aw,:as,:a,:gamNeg,:shape])
		al           = maxp[:,1:3]
		gam          = maxp[:,4:end]
	end

	@rput al
	@rput title
	@rput posteriors
	@rput analysis_folder

	R"""al = as.data.table(al)
			lbls = if(ncol(al) > 1){c(expression(paste('Posterior ',alpha[w])), expression(paste('Posterior ',alpha[s])),expression(paste('Posterior ',alpha)))}else{c(expression(paste('Posterior ',alpha)))}
			clrs = if(ncol(al) > 1){c('#30504f', '#e2bd9a', '#ab2710')}else{c('#ab2710')}

			if(nrow(al) == 1){
				tmp = as.data.table(posteriors[1])
				dal = suppressWarnings(data.table::melt(tmp[,1:3]))
				pal = ggplot(dal) + geom_density(aes(x=value,fill=variable),alpha=0.75) + scale_fill_manual('Posterior distribution',values=clrs ,labels=lbls) + theme_bw() + ggtitle(title) + xlab(expression(alpha)) + ylab("")
			}else{
				dal = suppressWarnings(data.table::melt(al))
				pal = suppressWarnings(ggplot(dal) + geom_density(aes(x=value,fill=variable),alpha=0.75) + scale_fill_manual('Posterior distribution',values=clrs ,labels=lbls) + theme_bw() + ggtitle(title) + xlab(expression(alpha)) + ylab(""))
			}
			suppressWarnings(ggsave(pal,filename=paste0(analysis_folder,'/map.png'),dpi=600))
			"""

	write(analysis_folder * "/map.tsv",maxp,delim='\t',header=true)

	return(posteriors,maxp)
end
