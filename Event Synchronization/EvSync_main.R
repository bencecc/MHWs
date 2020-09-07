#### ---- Main analysis of event synchronization. Significant pairwise ES between -------- ####
#### ---- any two sites/pixels is determined by comaprison with null model thresholds ---- ####
#### ---- given the same number of events (computed through EvSync_null_model). ---------- ####
#### ---- The function is called by Evsync_pipeline. ------------------------------------- ####
#### ------------------------------------------------------------------------------------- ####
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=36)
require(Rcpp)

# source cpp function to compute ES in parallel
sourceCpp('~/MHWs/code/EvSync_clust.cpp')
sourceCpp('~/MHWs/code/EvSync_sign_serial.cpp')

save.name <- 'es_sign'

args <- commandArgs(trailingOnly = TRUE)
dat <- mget(load(args[1]))[[1]]
tlen <- mget(load(args[2]))[[1]]
es_thr <- mget(load(args[3]))[[1]]
st <- as.numeric(args[4])
en <- as.numeric(args[5])

if(en==ncol(dat)) {
	en <- en-1
}

# Compute event synchronization bewteen selected columns through st:en and all other columns in the original data
system.time(es_mat <- EvSync_clust(dat[,st:en], dat[,st:ncol(dat)])) # the last column of es_mat is all zeros

# set vector of indices
vec <- st:en
IDs <- as.numeric(colnames(dat)) #preserve colnames (cell IDs) of input data

# select significant ESs; this id done by first identifying the number of MHW occurrences for each combination of sites/pixels:
# lx and ly; these become indices to select the threshold ES value derived from null models
system.time(es_sign_res <- foreach(i = 1:(ncol(es_mat)-1), .combine='rbind') %dopar% {
			
			tmp_res <- EvSync_sign_serial(i, es_mat[,i], as.matrix(es_thr), tlen, vec, IDs)
			tmp_res <- tmp_res[order(tmp_res[,1], tmp_res[,2]),]
			rows.rm <- which(rowSums(tmp_res)==0)
			tmp_res <- tmp_res[-rows.rm,]
			tmp_res
		
		})

if(dim(es_sign_res)[1] > 0) {
	
	colnames(es_sign_res) <- c("c1", "c2", "ESobs", "ESnull")
	
	assign(paste(save.name, st, '_', en, sep=""), value=es_sign_res, pos=1, inherits=T)
	outputName=paste(save.name, st, '_', en, ".RData",sep="")
	outputPath=file.path('~/EvSync_main', outputName)
	save(list=paste(save.name, st, '_', en, sep=""), file=outputPath)
	
}

assign(paste("es", st, '_', en, sep=""), value=es_mat, pos=1, inherits=T)
outputName=paste("es", st, '_', en, ".RData",sep="")
outputPath=file.path('~/EvSync_es', outputName)
save(list=paste("es", st, '_', en, sep=""), file=outputPath)

