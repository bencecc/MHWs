#### ---- Compute great circle distance from a coordinate file with cell number in the ---- ####
#### ---- first column and lon/lat in the second and thrid columns, respectively. --------- ####
#### ---- The function is called by Evsync_pipeline. -------------------------------------- ####
#### -------------------------------------------------------------------------------------- ####
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=72)
require(dplyr)

save.name <- 'gcd_'

args <- commandArgs(trailingOnly = TRUE)
coords <- mget(load(args[1]))[[1]]
st <- as.numeric(args[2])
en <- as.numeric(args[3])

if(en == nrow(coords)) {
	
	en = en - 1
}


system.time(dist_res <- foreach(i = st:en, .combine='rbind') %dopar% { 
			
			#cat('Doing iter ', i, ' of ', nrow(es_dat), '\n', sep = '')
			
			dist_dat <- foreach(j = (i+1):nrow(coords), .combine="rbind") %do% {
				
				gc_dist <- geosphere::distHaversine(c(coords[i,2], coords[i,3]), c(coords[j,2], coords[j,3]))
				
				cbind(c1=i, c2=j, dist=gc_dist/1000)
				
			}
			
			return(dist_dat)
		})

assign(paste(save.name, st, '_', en, sep=""), value=dist_res, pos=1, inherits=T)
outputName=paste(save.name, st, '_', en, ".RData",sep="")
outputPath=file.path('~/EvSync_all_dist', outputName)
save(list=paste(save.name, st, '_', en, sep=""), file=outputPath)


