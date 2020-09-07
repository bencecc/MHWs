#### ---- EVENT SYNCHRONIZATION PIPELINE ---- ####

# There are seven eight steps in the analysis:
# STEP 1: prepare a binary data frame of MHW occurrences for analysis: rows are days and columns are cells/pixels.
# STEP 2: compute the total number of events in each cell. The output is a vector.
# STEP 3: generate a two-column data frame of cell pairs with total number of MHW events in each cell;
#         each row of this data frame will be an entry for null model computation of event synchronization (ES)
#         given the total number of events in each of the two cells. The data frame includes unique combinations
#         of number of MHW events to avoid repeating the computation for a given pairwise combination of total
# 		  number of events more than once.
# STEP 4: run null models. The total number of events in each cell is randomized and ES is computed for 
#         each pair of cells 2000 times. A threshold value for ES is obtained at alpha=0.995 from null model
#         distributions; each threshold value will be the reference for observed ESs between cells with
#         the coresponding total number of MHW events in each cell. The output is a data frame of three columns;
# 		  the first two columns include the total number of MHW events between compared cells; the third column
#         is the threshold ES.
# STEP 5: identify significant links between cells. This is done by computing the observed ESs and select those
#         above threshold. The output is a matrix of four columns, the first 2 give the number of the cells for 
#         which ES has been computed, the third is the observed ES and the fourth is the threshold ES. This step
#         also saves the full matrix of ESs. 
# STEP 6: compute great-circle distance between cells with significant ES. The output is a vector of distances 
#         corresponding to the rows of the output matrix of step 5.
# STEP 7: probability density estimation
# STEP 8: compute great-circle distance between all selected cells in red_dat. The output is a matrix of distances. 
# STEP 9: remove folders.

# NOTE:   If required, taumax must be changed un cpp functions EvSync_nopar.cpp and EvSync_clust.cpp 
# IMPORTANT: Remember to make .sh scripts executable using chmod +x filename.sh.

#### ---- Help function to save files ---- ####
saveit <- function(..., string, file) {
	x <- list(...)
	names(x) <- string
	save(list=names(x), file=file, envir=list2env(x))
}

# -------------------------------------------- #

#### ---- Function to submit code in a termina editor as copy-paste and wait until a
#### ---- set of instructions is completed before proceeding to the next chunk of code  ---- ####
#### ---- NOTE: it works with PBS and the name of hpc must be provided: here is 'openhpc'    ---- ####
my_system_wait <- function(job_attr=paste('./mhw.sh', path_to_nc, sep=" "), ...) {
	
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

# --------------------------------------------------------------------------------------------- #

# require libraries
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=72)
require(dplyr)
require(data.table)
require(stringr)

# define folders for results 
setwd('~')
# main folder
if(!file.exists('EvSync_es'))
	dir.create('EvSync_es')
if(!file.exists('EvSync_tlen'))
	dir.create('EvSync_tlen')
if(!file.exists('EvSync_null'))
	dir.create('EvSync_null')
if(!file.exists('EvSync_main'))
	dir.create('EvSync_main')
if(!file.exists('EvSync_dist'))
	dir.create('EvSync_dist')
if(!file.exists('EvSync_rdist'))
	dir.create('EvSync_rdist')
if(!file.exists('EvSync_all_dist'))
	dir.create('EvSync_all_dist')

# ---- SET PATHS, FILES & DATES ---- #
outputPath <- "~/MHWs/data/"
output_name <- 'es_glob_1982_2018_'
setwd(outputPath)
tda_obj <- 'tda_glob.RData'
load(tda_obj)
red_dat <- tda_glob[, time:=NULL] # red_dat will become the focal dataset after declustering
hw_dat_name <- "hw_glob_df.RData" # hw_..._df contains coordinates and cell ids
hw_dat <- mget(load(paste(outputPath, "/", hw_dat_name, sep="")))[[1]]
time_period <- seq(as.Date("1982/1/1"), as.Date("2018/12/31"), by = "day")
#time_period <- seq(as.Date("1861/1/1"), as.Date("2005/12/31"), by = "day")
target_yrs <- as.character(1982:2018)

# --------------------------------------------------------------------------------------------- #

# prepare data for analysis: define time windows for event synchronization
# and subset the original file by the desired time window
tmp <- substr(time_period, 1,4)
row_sel <- which(tmp%in%target_yrs)
red_dat <- as.matrix(red_dat[row_sel, ])

#declustering
dec_fun <- function(x, ...) {
	d <- which(diff(x)==1)
	out <- double(length(x))
	out[d+1] <- 1
	return(out)
}

red_dat <- apply(red_dat, 2, dec_fun)
col.rm <- which(colSums(red_dat)<3) # select cells with at least three events
if(length(col.rm) > 0) {
		red_dat <- red_dat[,-col.rm]
	}

# ---- ARGUMENTS TO RUN FUNCTIONS ON CLUSTER ---- #
# for step 4
arg_tlen <- nrow(red_dat) #for analysis in step 4
arg_tlen_pairs <- file.path(outputPath, paste(output_name, 'tlen_pairs', ".RData", sep="")) 
# for step 5
arg_tlen_dat <- file.path(outputPath, paste(output_name, 'tlen', ".RData", sep="")) 
arg_dat <- file.path(outputPath, paste(output_name, 'dat', ".RData", sep=""))
arg_thr <- file.path(outputPath, paste(output_name, "thresholds", ".RData", sep=""))
# for step 6
arg_es_res <- file.path(outputPath, paste(output_name, "results", ".RData", sep=""))
arg_hw_dat <- file.path(outputPath, paste(output_name, "coords", ".RData", sep=""))

#### ---- Generate hw_df dataset containing only the cells present in red_dat
#### ---- for use in step 6 -------------------------------------------------- ####
hw_tmp <- hw_dat %>% filter(cell%in%as.numeric(colnames(red_dat)))
hw_tmp$x <- hw_tmp$x-179.5
hw_tmp <- hw_tmp %>% select(cell, x, y) %>% distinct()

# save coordinates of cells matching those selected in red_dat
saveit(dat=hw_tmp, string=paste(output_name, 'coords', sep=""),
		file=arg_hw_dat)
rm(hw_tmp)
rm(hw_dat)

#### ---------------------------------------------------------------------------------- ####

#### ----- STEP 1: SAVE THE TIME-SUBSETTED DATA ---- ####

saveit(dat=red_dat, string=paste(output_name, 'dat', sep=""),
		file=file.path(outputPath, paste(output_name, "dat", ".RData", sep="")))

#### ---- STEP 2: GET THE TOTAL NUMBER OF MHW EVENTS FOR EACH CELL ---- ####

tmp_events <- foreach(i = 1:ncol(red_dat), .combine='c') %do% {
	length(which(red_dat[,i]==1))
}

saveit(dat=tmp_events, string=paste(output_name, 'tlen', sep=""),
		file=file.path(outputPath, paste(output_name, "tlen", ".RData", sep="")))

#### ---- STEP 3: GENERATE UNIQUE PAIRWISE COMBINATIONS OF CELLS BASED ON TOTAL     ----- ####
#### ---- NUMBER OF MHW EVENTS; THIS IS DONE TO MINIMIZE THE NUMBER OF NULL MODELS  ----- ####

# generate indices to partition analysis among cores
ln1 <- seq(1, length(tmp_events), by=1000)
ln2 <- seq(1000, length(tmp_events), by=1000)
evsync_tlen_parms.txt  <- cbind(ln1[1:length(ln1)-1], c(ln2[1:length(ln2)-1], length(tmp_events)))
write.table(evsync_tlen_parms.txt, file='~/MHWs/data/evsync_tlen_parms.txt', sep="\t", row.names=F, col.names=F)

# run in parallel
my_system_wait(job_attr=paste('~/MHWs/code/./evsync_tlen.sh', arg_tlen_dat, sep=" "))

setwd('~/EvSync_tlen')
file_list <- as.list(list.files())
es_pairs_tmp <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	mget(load(file_list[[i]]))[[1]]
}

es_pairs <- es_pairs_tmp %>% distinct(c1,c2, .keep_all=T)
saveit(dat=es_pairs, string=paste(output_name, 'tlen_pairs', sep=""),
		file=file.path(outputPath, paste(output_name, "tlen_pairs", ".RData", sep="")))

#### ---- STEP 4: RUN NULL MODELS TO GENERATE THRESHOLDS VALUES OF ES GIVEN ----- ####
#### ----       THE TOTAL NUMBER OF MHW OCCURRENCES BETWEEN ANY TWO CELLS           ----- ####
ln1 <- seq(1, nrow(es_pairs), by=500)
ln2 <- seq(500, nrow(es_pairs), by=500)
evsync_null_parms.txt  <- cbind(ln1[1:length(ln1)-1], c(ln2[1:length(ln2)-1], nrow(es_pairs)))
write.table(evsync_null_parms.txt, file='~/MHWs/data/evsync_null_parms.txt', sep="\t", row.names=F, col.names=F)

my_system_wait(job_attr=paste('~/MHWs/code/./evsync_null.sh', arg_tlen_pairs, arg_tlen, sep=" "))

setwd('~/EvSync_null')
file_list <- as.list(list.files())
es_null <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	mget(load(file_list[[i]]))[[1]]
}

saveit(dat=es_null, string=paste(output_name, 'thresholds', sep=""),
		file=file.path(outputPath, paste(output_name, "thresholds", ".RData", sep="")))

#### ---- STEP 5: IDENTIFY SIGNIFICANT LINKS BETWEEN CELLS ----- ####
ln1 <- seq(1, length(tmp_events), by=500)
ln2 <- seq(500, length(tmp_events), by=500)
evsync_main_parms.txt  <- cbind(ln1[1:length(ln1)-1], c(ln2[1:length(ln2)-1], length(tmp_events)))
write.table(evsync_main_parms.txt, file='~/MHWs/data/evsync_main_parms.txt', sep="\t", row.names=F, col.names=F)

my_system_wait(job_attr=paste('~/MHWs/code/./evsync_main.sh', arg_dat, arg_tlen_dat, arg_thr, sep=" "))

setwd('~/EvSync_main')
file_list <- as.list(list.files())
es_sign <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	mget(load(file_list[[i]]))[[1]]
}

es_sign <- es_sign[order(es_sign[,"c1"], es_sign[,"c2"]), ]

saveit(dat=es_sign, string=paste(output_name, 'results', sep=""),
		file=file.path(outputPath, paste(output_name, "results", ".RData", sep="")))

#### ---- Save full Event Synchronization matrix ---- ####
setwd('~/EvSync_es')

# ensure that the first file in file_list is es1_500
file_list <- as.list(list.files())
m1 <- mget(load(file_list[[1]]))[[1]]
rr <- nrow(m1)
cc <- ncol(m1)
m <- matrix(0, nrow=rr, ncol=rr)

# order file_list
ln1 <- seq(1, rr, by=800)
ln2 <- seq(800, rr, by=800)
mysort <- paste('es', ln1[1:length(ln1)-1], "_", c(ln2[1:length(ln2)-1], rr), '.RData', sep="") #order is important
file_list <- file_list[order(charmatch(file_list, mysort))]

for(i in 1:length(file_list)) {
	if(i == 1) {
		m[i:rr, 1:cc] <- m1
	}
	else {
		tmp_dat <-  mget(load(file_list[[i]]))[[1]]
		delta_row <- rr-nrow(tmp_dat)
		m[(delta_row+1):rr, ((i-1)*cc+1):((i-1)*cc+ncol(tmp_dat))] <- tmp_dat
	}
}

saveit(dat=m, string=paste(output_name, 'es', sep=""),
		file=file.path(outputPath, paste(output_name, "es", ".RData", sep="")))

#### ---- STEP 6: COMPUTE great-circle DISTANCE BETWEEN CELLS INVOLVED IN SIGNIFICANT LINKS ---- ####
ln1 <- seq(1, nrow(es_sign), by=200000)
ln2 <- seq(200000, nrow(es_sign), by=200000)
evsync_dist_parms.txt  <- cbind(ln1[1:length(ln1)-1], c(ln2[1:length(ln2)-1], nrow(es_sign)))
write.table(evsync_dist_parms.txt, file='~/MHWs/data/evsync_dist_parms.txt', sep="\t", row.names=F, col.names=F)

my_system_wait(job_attr=paste('~/MHWs/code/./evsync_dist.sh', arg_es_res, arg_hw_dat, sep=" "))

setwd('~/EvSync_rdist')
file_list <- as.list(list.files())
es_dist <- foreach(i=1:length(file_list), .combine='rbind') %dopar% {
	mget(load(file_list[[i]]))[[1]]
}

es_dist <- es_dist[order(es_dist[,"c1"], es_dist[,"c2"]),]
		
saveit(dat=es_dist, string=paste(output_name, 'rdist', sep=""),
		file=file.path(outputPath, paste(output_name, "rdist", ".RData", sep="")))

#### ---- STEP 7: Probability Density Estimation ---- ####
require(kdensity)

sig_dist_mat <- es_dist
set.seed(12345)
sig_dist_samp <- sig_dist_mat[sample(1:nrow(sig_dist_mat), 1e6), "dist"]
sig_dist_dens <- kdensity(sig_dist_samp, na.rm=T,
		start = "gumbel", kernel = "gaussian")
space <- seq(50, 20000, length.out=2000)
sig_gcdist_density=sig_dist_dens(space)

saveit(dat=sig_gcdist_density, string=paste("es","_","declust","_",target_yrs[1],"_",target_yrs[length(target_yrs)], sep=""),
		file=file.path(outputPath, paste("es","_","declust","_",target_yrs[1],"_",target_yrs[length(target_yrs)], ".RData", sep="")))

#### ---- STEP 8: COMPUTE great-circle DISTANCE BETWEEN ALL CELLS     ---- ####
#### ---- To do only once with red_dat set as the full tda object;    ---- ####
#### ---- thus, repeat steps 1-7 omitting code from lines 111 to 127. ---- ####
#### ---- Step 7 will provide the KDE for all possible GCDs, while    ---- ####
#### ---- the code below will save the full distance matrix.          ---- ####

# uses the same parms indices as in step 5: evsync_main_parms.txt
my_system_wait(job_attr=paste('~/MHWs/code/./evsync_all_dist.sh', arg_hw_dat, sep=" "))

setwd('~/EvSync_all_dist')

# order files in file_list starting from gcd_1_500
file_list <- as.list(list.files())
# m1 is defined in step 5
rr <- nrow(m1)
d <- matrix(0, nrow=rr, ncol=rr)

# order file_list
ln1 <- seq(1, rr, by=500)
ln2 <- seq(500, rr, by=500)
mysort <- paste('gcd_', ln1[1:length(ln1)-1], "_", c(ln2[1:length(ln2)-1], rr), '.RData', sep="") #order is important
file_list <- file_list[order(charmatch(file_list, mysort))]

for(i in 1:length(file_list)) {
	tmp_dat <-  mget(load(file_list[[i]]))[[1]]
	for(j in 1:nrow(tmp_dat)) {
		d[tmp_dat[j,2], tmp_dat[j,1]] <- tmp_dat[j,3]
	}
}

saveit(dat=d, string=paste(output_name, 'all_dist', sep=""),
		file=file.path(outputPath, paste(output_name, "all_dist", ".RData", sep="")))

#### ---- STEP 9: remove temporary folders ---- ####
setwd('~')
if(file.exists('EvSync_es'))
	unlink('EvSync_es', recursive=T)
if(file.exists('EvSync_tlen'))
	unlink('EvSync_tlen', recursive=T)
if(file.exists('EvSync_null'))
	unlink('EvSync_null', recursive=T)
if(file.exists('EvSync_main'))
	unlink('EvSync_main', recursive=T)
if(file.exists('EvSync_dist'))
	unlink('EvSync_dist', recursive=T)
if(file.exists('EvSync_rdist'))
	unlink('EvSync_rdist', recursive=T)
if(file.exists('EvSync_all_dist'))
	unlink('EvSync_all_dist', recursive=T)

#### ---- END OF ANALYSIS ---- ####

