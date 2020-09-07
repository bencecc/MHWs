#### ---- Run null models of event synchronization of number of MHW occurrences ------- ####
#### ---- obtained from function EvSync_tlen. These are pairwise number of events ----- ####
#### ---- for each combination of sites/pixels for which ES need to be computed. ------ ####
#### ---- The output of EvSync_null_model is the ES threshold for each combination ---- ####
#### ---- of occurrences to identify significance links betwen sites/pixels in the ---- ####
#### ---- analysis.The function is called by Evsync_pipeline. ------------------------- ####
#### ---------------------------------------------------------------------------------- ####
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=72)
require(Rcpp)

# source cpp function to compute ES from two vectors of indices of ES occurrences
sourceCpp('~/MHWs/code/EvSync_nopar.cpp')

save.name <- 'es_q'
nrand <- 2000

args <- commandArgs(trailingOnly = TRUE)
dat <- mget(load(args[1]))[[1]]
tlen <- as.numeric(args[2])
st <- as.numeric(args[3])
en <- as.numeric(args[4])

q_res <- foreach(i = st:en, .combine='rbind') %dopar% {
	
	lx=dat[i,1]
	ly=dat[i,2]
	
	es <- foreach(k = 1:nrand, .combine='c') %do% {
		
		ex <- rep(0, tlen)
		ey <- rep(0, tlen)
		ex[sample(1:tlen, lx)] <- 1
		ey[sample(1:tlen, ly)] <- 1
		ind.ex <- which(ex==1, arr.ind=T)
		ind.ey <- which(ey==1, arr.ind=T)
		
		EvSync_nopar(ind.ex,ind.ey)
		
	}
	
	q <- as.numeric(quantile(es, 0.995))
	as.matrix(data.frame(c1=lx, c2=ly, q=q))
}

assign(paste(save.name, st, '_', en, sep=""), value=q_res, pos=1, inherits=T)
outputName=paste(save.name, st, '_', en, ".RData",sep="")
outputPath=file.path('~/EvSync_null', outputName)
save(list=paste(save.name, st, '_', en, sep=""), file=outputPath)


