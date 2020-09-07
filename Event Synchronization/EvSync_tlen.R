#### ---- Find unique combinations of number of MHW occurrences between pairwise ---- ####
#### ---- cells/pixels from a tda object to provide input to null models ------------ #### 
#### ---- in event synchronization analysis. The function is called by -------------- ####
#### ---- Evsync_pipeline. ---------------------------------------------------------- ####
#### -------------------------------------------------------------------------------- ####
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=72)
require(dplyr)

save.name <- 'tlen'

args <- commandArgs(trailingOnly = TRUE)
dat <- mget(load(args[1]))[[1]]
st <- as.numeric(args[2])
en <- as.numeric(args[3])

if(en==length(dat))
	en = en - 1

ev_pairs_tmp <- foreach(i = st:en, .combine='rbind') %dopar% {
	
	out <- foreach(j = (i+1):length(dat), .combine='rbind') %do% {
		
		cbind(dat[i], dat[j])
		
	}
	
	colnames(out) <- c("c1","c2")
	
	as.data.frame(out) %>% distinct(c1,c2, .keep_all=T)
}

ev_pairs <- ev_pairs_tmp %>% distinct(c1, c2, .keep_all=T)

assign(paste(save.name, st, '_', en, sep=""), value=ev_pairs, pos=1, inherits=T)
outputName=paste(save.name, st, '_', en, ".RData",sep="")
outputPath=file.path('~/EvSync_tlen', outputName)
save(list=paste(save.name, st, '_', en, sep=""), file=outputPath)
