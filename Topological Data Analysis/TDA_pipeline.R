#### ---- PIPLINE FOR A FULL TOPOLOGICAL DATA ANALYSIS OF MHWs ---- ####
#### ---- The script uses a .nc file with remotely sensed or ---- ####
#### ---- simulated SST timeseries for the global ocean and ------- ####
#### ---- runs a full topological data analysis on the data.------- ####
#### ---- Alternatively, one can start the analysis from Step ----- ####
#### ---- 3 using the provided tda and hw files. ------------------ ####

#### ---- IMPORTANT: Remember to make .sh scripts executable ------ ####
#### ----- using chmod +x filename.sh. ---------------------------- ####

#### ---- Names for input/output files and output folders---- ####
# file identifier
fi_id <- 0
# name for data tibble containing all MHW info and precursor for covariate file
hw_output_name <- 'hw_glob_'
# name for tda output file
tda_name <- 'tda_glob_'

# define folders for results 
setwd('~')
# main folder
if(!file.exists('tda_results'))
	dir.create('tda_results')
setwd('~/tda_results')
# specific for individual nc file
if(!file.exists(paste(tda_name, fi_id, '_out', sep="")))
	dir.create(paste(tda_name, fi_id, '_out', sep=""))
# output path
outputPath=file.path(paste('~/tda_results/', paste(tda_name, fi_id, '_out', sep=""), sep=""))
# set path to nc file
path_to_nc <- file.path(paste('~/MHWs/data/', 'sst_', fi_id, '.nc', sep=""))

#### ---- Require R libraries ---- ####
require(TDAmapper)
require(igraph)
require(reticulate)
require(dplyr)
require(tidyr)
require(abind)
require(data.table)
require(stringr)

require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=30) #adjust for the resources at your disposal

#### ------------------------------------------------------------------------ ####
#### ---- Function to save files -------------------------------------------- ####
saveit <- function(..., string, file) {
	x <- list(...)
	names(x) <- string
	save(list=names(x), file=file, envir=list2env(x))
}

#### ---- Function to wait until code is completed on hpc to continue wit ---- ####
#### ---- the function. ------------------------------------------------------ ####
#### ---- NOTE: it works with PBS and name of hpc needs to be set: ----------- ####
#### ---- here is 'openhpc'. ------------------------------------------------- ####
my_system_wait <- function(job_attr=paste('~/MHWs/code/./mhw.sh', path_to_nc, sep=" "), ...) {
	
	sys_out <- system(job_attr, wait=TRUE, intern=TRUE)
	n_collect <- list()
	for(k in 1:length(sys_out)) {
		n_tmp <- charmatch('JobID = ', sys_out[k])
		if(!is.na(n_tmp)) {
			str.tmp <- strsplit(sys_out[k], split=c('JobID','.'))[[1]][2]
			n_collect[[k]] <- substr(str.tmp, start=regexpr('=', str.tmp)[1]+2,
					stop=regexpr('.openhpc', str.tmp)[1]-1)
			
		}
	}
	
	n_collect <- unlist(plyr::compact(n_collect))
	
	q_out <- system('qstat -a', wait=TRUE, intern=T)
	running_jobs <- sapply(n_collect, function(x) str_extract(q_out, x))
	while(length(dimnames(running_jobs)[[2]])>0) {
		Sys.sleep(60)
		q_out <- system('qstat -a', wait=TRUE, intern=T)
		running_jobs <- sapply(n_collect, function(x) str_extract(q_out, x))
	} 
}

#### ---- STEP 1: detect MHW from .nc files; multicore computing using function ----------- ####
#### ---- HeatWaves.R HeatWavesClim.R - the latter uses climatological -------------------- ####
#### ---- historical data (1861-2005) and is employed with rcp26, rcp85 and, -------------- ####
#### ---- optionally, with contemporary SST satellite data. Function name, ---------------- ####
#### ---- HeatWaves.R or HeatWavesClim.R, must be set in mhw.sh; also set ----------------- ####
#### ---- focal_var as'tos' when using historical/rcp85/rcp26 .nc files or as ------------- ####
#### ---- 'SST' with remotely sensed SSTs. If tos is used, temperature unit is ------------ ####
#### ---- in kelvin and one may want to subtract 273.15 at line 71 of HeatWaves.R --------- ####
#### ---- or at line 90 of HeatWavesClim.R. Also comment/uncomment lines 110-153 ---------- ####
#### ---- in HeatWavesClim.R depending on whether the analysis uses simulated ------------- ####
#### ---- (rcp85/rcp26) or remotely sensed SST data. Script HeatWaves.R is ---------------- ####
#### ---- currently set for an analysis of contemporary satellite MHWs (sst variable). ---- ####

setwd('~')
if(!file.exists('tmp_files'))
	dir.create('tmp_files')

mhw_parms_1 <- 1:180 #separate analysis by latitude
write.table(mhw_parms_1, file='~/MHWs/data/mhw_parms_1.txt', sep="\t", row.names=F, col.names=F)
my_system_wait(job_attr=paste('~/MHWs/code/./mhw.sh', path_to_nc, sep=" "))

#### ---- STEP 2: generate TDA and initial covariate files; a TDA file has time ---- ####
#### ---- as the first column and geographic pixels as the other columns; data  ---- ####
#### ---- entries are occurrences (0,1) of MHWs. ----------------------------------- ####   
setwd('~/tmp_files')

file_list <- as.list(list.files())

mysort <- paste('lat_', 1:length(file_list), '.RData', sep="") #if order is important
file_list <- file_list[order(charmatch(file_list, mysort))]

file_df <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	mget(load(file_list[[i]]))[[1]]
}

# save heatwave file
saveit(dat=file_df, string=paste(hw_output_name, fi_id, '_df', sep=""),
		file=file.path(outputPath, paste(hw_output_name, fi_id, "_df", ".RData", sep="")))

# generate tda file
cells <- unique(file_df$cell)
tmp_df <- foreach (i=1:length(cells), .combine = 'rbind') %dopar% {
	cell_df <- file_df %>% filter(cell%in%cells[i])
	
	out <- apply(cell_df, 1, function(x) {
				nrep <- x['duration']
				sub_time <- x['index_start']:x['index_end']
				sub_df <- data.frame(cell=as.numeric(rep(x['cell'], nrep)),
						time=sub_time, #int_cum=as.numeric(rep(x['intensity_cumulative'], nrep)),
						hw_pa=as.numeric(rep(1, nrep)))
			})
	
	bind_rows(out)
}

tda_df <- as.data.table(tmp_df)
rm(tmp_df)

# split analysis in two files to reduce size
target.cell <- tda_df[ceiling(nrow(tda_df)/2),'cell']
cel.sel <- which(tda_df$cell%in%target.cell)

tmp <- tda_df[1:cel.sel[length(cel.sel)],]
tmp1 <- tda_df[(cel.sel[length(cel.sel)]+1):nrow(tda_df),]

tda_df11 <- dcast(tmp, time~cell, value.var='hw_pa', fill=0)
tda_df12 <- dcast(tmp1, time~cell, value.var='hw_pa', fill=0)
tda_final <- merge(tda_df11, tda_df12, by='time')
#save tda file
saveit(dat=tda_final, string=paste(tda_name, fi_id, sep=""),
		file=file.path(outputPath, paste(tda_name, fi_id, ".RData", sep="")))

rm(tda_df11)
rm(tda_df12)
rm(tmp)
rm(tmp1)
rm(target.cell)
rm(cells)
rm(tda_df)
gc()

#### ---- STEP 3: generate covariates; tda_final is one of the TDA objects   ------------ ####
#### ---- (global, historical, rcp85, rcp26) containing data on presence/absence -------- ####
#### ---- of MHWs resulting from the analysis above. These files are  provided
#### ---- in the repository.
#### ------------------------------------------------------------------------------------ ####

# generate indicator vars to split the analysis across multiple cores;
len1 <- ceiling(seq(1, nrow(tda_final), length.out=7))
len2 <- len1[2:6]
parms_cov_1 <- cbind(c(len1[1],len1[2:6]+1), c(len2, len1[7]))
write.table(parms_cov_1, file='~/MHWs/data/parms_cov_1.txt', sep="\t", row.names=F, col.names=F)

setwd('~')
if(!file.exists('tmp_cov_res'))
	dir.create('tmp_cov_res')

arg_cov1 <- file.path(outputPath, paste(tda_name, fi_id, ".RData", sep=""))
arg_cov2 <- file.path(outputPath, paste(hw_output_name, fi_id, "_df", ".RData", sep=""))

my_system_wait(job_attr=paste('~/MHWs/code/./covar_gen.sh', arg_cov1, arg_cov2, sep=" "))

setwd('~/tmp_cov_res')
file_list <- as.list(list.files())
#order
catch_files <- paste('cov_', parms_cov_1[,1], '.RData', sep="")
file_list <- file_list[order(charmatch(file_list, catch_files))]
tda_covar <- foreach(i=1:length(file_list), .combine=rbind) %do% {
	mget(load(file_list[[i]]))[[1]]
	
}

# save covar file
saveit(dat=tda_covar, string=paste(tda_name, fi_id, "_covar", sep=""),
		file=file.path(outputPath, paste(tda_name, fi_id, "_covar", ".RData", sep="")))

#### ---- STEP 4: optimization of tda analysis on parameter grid (res, gain) ---- ####
#### ---- The resulting tda_res object is a list that includes the adjacency ---- ####
#### ---- matrix used to build networks and to compute network statistics.   ---- ####
require(uwot)
require(bigstatsr)
# generate lens: umap approach
fmb_rand_dat <- as_FBM(tda_final)
red_dat <- big_randomSVD(fmb_rand_dat, k=50, ncores=30)
pca_scores <- predict(red_dat)
umap_res <- umap(pca_scores)
rm(red_dat)
rm(fmb_rand_dat)
saveit(dat=umap_res, string=paste(tda_name, fi_id, "_umap", sep=""),
		file=file.path(outputPath, paste(tda_name, fi_id, "_umap", ".RData", sep="")))

# run optimization:
setwd('~')
if(!file.exists('tmp_optim'))
	dir.create('tmp_optim')

arg_tda <- file.path(outputPath, paste(tda_name, fi_id, ".RData", sep=""))
arg_umap <- file.path(outputPath, paste(tda_name, fi_id, '_umap', ".RData", sep=""))
arg_covar <- file.path(outputPath, paste(tda_name, fi_id, '_covar', ".RData", sep=""))

my_system_wait(job_attr=paste('~/MHWs/code/./tda_optim.sh', arg_tda, arg_umap, arg_covar, sep=" "))

setwd('~/tmp_optim')
file_list <- as.list(list.files())
#Extract Fstat, res and gain from data files
wrap_dat <- abind(lapply(sapply(file_list, function(x) 
						{tmp_dat <- mget(load(x))
							tmp_dat},
						simplify = TRUE), function(x) {
					tda_tmp <- x[[2]]
					naratio <- length(which(unlist(tda_tmp$points_in_level)==0))/length(tda_tmp$points_in_level) #ratio of empty nodes
					singratio <- length(which(unlist(lapply(tda_tmp$points_in_vertex, function(x) length(x)))==1))/length(tda_tmp$points_in_vertex) #ratio of singleton vertices
					dim.adj <- dim(tda_tmp$adjacency)[1]
					c(x[[1]], dim.adj=dim.adj, NA_Ratio=naratio, SingRatio=singratio)}), along=0) 

# identify the tda_output that maximises differences in cumulative MHW
# intensity among nodes compared to variability within nodes
wrap_df <- data.frame(wrap_dat)
wrap_df$data.name <- rownames(wrap_dat)
optim_filter <- data.frame(wrap_df)
optim_df <- optim_filter %>% mutate(balance=1/scale(Fstat)) %>%
		slice(which.min(balance))
sel <- which(unlist(file_list)%in%paste(optim_df['data.name'], '.RData', sep=""))

tda_res <- mget(load(unlist(file_list)[sel]))[[1]][[2]]
tda_res <- append(tda_res, optim_df)
saveit(dat=tda_res, string=paste(tda_name, fi_id, "_res", sep=""),
		file=file.path(outputPath, paste(tda_name, fi_id, "_res", ".RData", sep="")))

#### ---- STEP 5: compute network statistics, temporal connectivity matrix and degree ------------------------------- ####
#### ---- use function temporal_connectivity to obtain a list where the first element is the ------------------------ ####
#### ---- temporal connectivity matrix and the second element is the nomalized temporal degreee vector. ------------- ####
#### ---- IMPORTANT: module cpalgorithm is not compatible with conda - you may need to force reticulate ------------- ####
#### ---- to look for a python version outside conda using Sys.setenv(RETICULATE_PYTHON = "your_path_to_python") ---- ####
#### ---- and install using py_install("cpalgorithm"). -------------------------------------------------------------- ####
source('~/MHWs/code/net_stats.R')
source('~/MHWs/code/temporal_degree.R')
# import pyton modules for core-perifery network analysis
cp <- import('cpalgorithm', as='cp', delay_load=F)
nx <- import('networkx', as='nx', delay_load=F)

tda_stats <- net_stats(tda_res)
tcm_deg <- temporal_degree(tda_res)
#normalize degree
tcm_deg <- (tcm_deg-min(tcm_deg))/(max(tcm_deg)-min(tcm_deg))

saveit(dat=tda_stats, string=paste(tda_name, fi_id, "_stats", sep=""),
		file=file.path(outputPath, paste(tda_name, fi_id, "_stats", ".RData", sep="")))
saveit(dat=tcm_deg, string=paste(tda_name, fi_id, "_degree", sep=""),
		file=file.path(outputPath, paste(tda_name, fi_id, "_degree", ".RData", sep="")))

#### ---- STEP 6: run null models ----------------------------------------------------------- ####
#### ---- NOTE: This analysis is very slow: 1000 simulations require one week of ------------ ####
#### ---- computation on an hpc of eight nodes and 36 cores for each node. ------------------ #### 
#### ---- NOTE: to modify the number of simulations change file nul_net_$ in run_null.sh ---- #### 
setwd('~')
if(!file.exists('tmp_null_res'))
	dir.create('tmp_null_res')

null_net_1 <- 1:1000 #random iterations
write.table(null_net_1, file='~/MHWs/data/null_net_1.txt', sep="\t", row.names=F, col.names=F)

arg_res <- file.path(outputPath, paste(tda_name, fi_id, '_res', ".RData", sep=""))
my_system_wait(job_attr=paste('~/MHWs/code/./run_null.sh', arg_tda, arg_res, sep=" "))

setwd('~/tmp_null_res')
file_list <- as.list(list.files())

obs <- NULL
null_phase <- NULL
null_rewire <- NULL
null_degree <- list()

for (i in 1:length(file_list)) {
	dat <- mget(load(file_list[[i]]))[[1]]
	obs <- rbind(obs, dat[[1]])
	null_phase <- rbind(null_phase, dat[[2]])
	null_rewire <- rbind(null_rewire, dat[[3]])
	null_degree[[i]] <- as.data.frame(matrix(dat[[4]], nrow=1))	

}

require(data.table)

null_degree <- data.table::rbindlist(null_degree, fill=T)

tda_null_stats <- list(obs, null_phase, null_rewire, null_degree)
saveit(dat=tda_null_stats, string=paste(tda_name, fi_id, "_null_stats", sep=""),
		file=file.path(outputPath, paste(tda_name, fi_id, "_null_stats", ".RData", sep="")))

#### ---- STEP 7: remove temporary folders ---- ####
setwd('~')
if(file.exists('tmp_files'))
	unlink('tmp_files', recursive=T)
if(file.exists('tmp_cov_res'))
	unlink('tmp_cov_res', recursive=T)
if(file.exists('tmp_optim'))
	unlink('tmp_optim', recursive=T)
if(file.exists('tmp_null_res'))
	unlink('tmp_null_res', recursive=T)

#### ---- END OF ANALYSIS ---- ####

