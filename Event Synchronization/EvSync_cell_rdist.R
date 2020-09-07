#### ---- Obtain distances between significantly synchronized sites. ---- ####
#### ---- The function is called by Evsync_pipeline. -------------------- ####
#### -------------------------------------------------------------------- ####
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=18)
require(dplyr)

save.name <- 'gcd_'

args <- commandArgs(trailingOnly = TRUE)
es_dat <- mget(load(args[1]))[[1]]
hw_dat <- mget(load(args[2]))[[1]]
st <- as.numeric(args[3])
en <- as.numeric(args[4])

dist_res <- foreach(i = st:en, .combine='rbind') %dopar% { 
	
	coords=hw_dat %>% filter(cell%in%es_dat[i,'c1']|cell%in%es_dat[i,'c2'])
	gc_dist <- geosphere::distHaversine(c(coords[1,2], coords[1,3]), c(coords[2,2], coords[2,3])) 
	cbind(c1=round(coords[1,1], 0), c2=round(coords[2,1], 0), dist=gc_dist/1000)
	
}


assign(paste(save.name, st, '_', en, sep=""), value=dist_res, pos=1, inherits=T)
outputName=paste(save.name, st, '_', en, ".RData",sep="")
outputPath=file.path('~/EvSync_rdist', outputName)
save(list=paste(save.name, st, '_', en, sep=""), file=outputPath)


