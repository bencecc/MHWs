# Data driven parameter search to select optimal resolution (intervals)
# and gain (percent overlap) parameters in TDA analysis. The function is called
# by TDA_pipeline.R
###############################################################################

source('~/MHWs/code/mapper2D_parallel.R')

require(TDAmapper)
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=20)

args <- commandArgs(trailingOnly = TRUE)

my_dat <- mget(load(args[1]))[[1]]
my_dat <- my_dat[,2:ncol(my_dat)]
test_umap <- mget(load(args[2]))[[1]]
covar_dat <- mget(load(args[3]))[[1]]

res <- as.numeric(args[4])
gain <- as.numeric(args[5])

## save name
save.name <- 'optim_'

###############################################################################
system.time(out_mapper <- mapper2D_parallel(
				data=my_dat,
				is.distance=F,
				ncores=1,
				filter_values=list(test_umap[,1], test_umap[,2]),
				num_intervals = c(res,res),
				percent_overlap=gain,
				num_bins_when_clustering = 10))

nodes <- mapperVertices(out_mapper, 1:nrow(my_dat))

#nodes_int <- apply(nodes, 1, function(x) {
system.time(nodes_int <- foreach(i=1:nrow(nodes)) %dopar% {
			x <- nodes[i,]
			#identify the rows of the time_dat contributing to the node
			str.tmp <- strsplit(as.character(x[,'Nodename']), split=c('V',':'))[[1]][2]
			substr(str.tmp, start=regexpr(':', str.tmp)[1], stop=regexpr(':', str.tmp)[1]) <- ','
			focal_times=as.vector(as.numeric(unlist(strsplit(str.tmp,","))))[-1]
			sub_dat <- covar_dat[focal_times,'intensity_cumulative']
			list(sub_dat)
		})

score_tmp <- sapply(nodes_int, function(x) {
			c(mean(x[[1]]), var(x[[1]]))
		})

score <- var(score_tmp[1,], na.rm=T)/mean(score_tmp[2,], na.rm=T)

out <- list(c(res=res, gain=gain, Fstat=score), out_mapper=out_mapper)

# save output
assign(paste(save.name, res, '_', gain, sep=""), value=out, pos=1, inherits=T)
outputName=paste(save.name, res, '_', gain, ".RData",sep="")
outputPath=file.path('~/tmp_optim', outputName)
save(list=paste(save.name, res, '_', gain, sep=""), file=outputPath)
