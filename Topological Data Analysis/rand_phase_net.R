# Null model randomization of network and TDA analysis. The function is called
# by TDA_pipeline.R

###############################################################################
source('~/MHWs/code/mapper2D_parallel.R')
source('~/MHWs/code/net_stats.R')
source('~/MHWs/code/temporal_degree.R')

require(uwot)
require(TDAmapper)
require(fastcluster)
require(Matrix)
require(wordspace)
require(VertexSort)
require(data.table)
require(bigstatsr)
require(igraph) 
require(fractal)
require(reticulate)

cp <- import('cpalgorithm', as='cp', delay_load=F)
nx <- import('networkx', as='nx', delay_load=F)

args <- commandArgs(trailingOnly = TRUE)
dat <- mget(load(args[1]))[[1]]
dat <- dat[,2:ncol(dat)]
dat_res <- mget(load(args[2]))[[1]]
iter <- as.numeric(args[3])

res <- dat_res$res
gain <- dat_res$gain

i_graph <- graph.adjacency(dat_res$adjacency, mode="undirected")

obs_stats <- net_stats(dat_res)
obs_degree <- temporal_degree(dat_res) #this version of function temporal_deree returns the unnormalized degree;

require(foreach, quietly=T)
library(parallel)
cl <- makeCluster(36)
doParallel::registerDoParallel(cl)

# set integer for constant phase randomization in function surrogate; otherwise seed=0
#my.seed <- as.numeric(Sys.time())

system.time(rand_dat <- foreach(col=iterators::iter(dat, by="col"),
						.combine='c', .packages=('fractal')) %dopar% {
					num <- as.numeric(col)
					out=as.vector(surrogate(num, method="aaft", sdf=NULL, seed=0))
					list(out)
				})

stopCluster(cl)
rm(dat)
gc()

rand_dat <- matrix(unlist(rand_dat), ncol=length(rand_dat), byrow=F)

fmb_rand_dat <- as_FBM(rand_dat)
red_dat <- big_randomSVD(fmb_rand_dat, k=50, ncores=30)
pca_scores <- predict(red_dat)
test_umap <- umap(pca_scores)
rm(red_dat)
rm(fmb_rand_dat)

rand_dat <- as.data.table(rand_dat)

system.time(out_mapper <- mapper2D_parallel(
				data=rand_dat,
				is.distance=F,
				ncores=1,
				filter_values=list(test_umap[,1], test_umap[,2]),
				num_intervals = c(res,res),
				percent_overlap=gain,
				num_bins_when_clustering = 10))

gc()

rand_phase <- net_stats(out_mapper)

rand_graph <- rewire(i_graph, keeping_degseq()) #dpr(i_graph, 1)
rand_rewire <- net_stats(rand_graph)

#degree of temporal connectivity matrix
# mapper objects generated under a constant phase null model are very large matrices and computing the 
# node degree for the temporal connectivity matrix of these objects requires a considerable amout of RAM
# to obivate this, one may rin the code omitting the two lines below for constamt phase null models
rand_degree <- temporal_degree(out_mapper) #this version of function temporal_deree returns the un-normalized degree;
norm_rand_degree <- (rand_degree-min(rand_degree))/(max(obs_degree)-min(rand_degree)) #normalization with respect to observed degree
out <- list(obs_stats=obs_stats, rand_phase=rand_phase, rand_rewire=rand_rewire, rand_degree=norm_rand_degree)

outputName=paste('r_', iter, ".RData", sep="")
outputPath=file.path('~/tmp_null_res', outputName)
assign(paste('r_', iter, sep=""), value=out, pos=1, inherits=T)
save(list=paste('r_', iter, sep=""), file=outputPath)


