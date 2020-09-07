# Generate covariates to color tda network graphs; function called
# by TDA_pipeline.R
###############################################################################
require(dplyr)
require(data.table)
require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=20)

args <- commandArgs(trailingOnly = TRUE)
tda_dat <- mget(load(args[1]))[[1]]
hw_dat_df <- mget(load(args[2]))[[1]]
time_first <- as.numeric(args[3])
time_last <- as.numeric(args[4])

time_sel <- unique(tda_dat[time_first:time_last, 'time'])

out <- foreach(i=1:nrow(time_sel), .combine=rbind) %dopar% {
	
	
	dat <- tda_dat %>% filter(time%in%time_sel[i])
	cel_sel <- dat %>% .[,colSums(.)>0] %>% select(., -time) %>% 
			colnames(.) %>% as.numeric(.)
	cross_time <- as.numeric(dat[,'time']) 
	sub_dat <- hw_dat_df %>% filter(cell%in%cel_sel)
	out_df <- sub_dat %>%
			do(data.frame(sel=data.table::between(cross_time, .$index_start, .$index_end))) %>%
			do(data.frame(sel=which(.$sel==T))) %>%
			do(data.frame(sub_dat[.$sel,])) %>%
			summarise(duration=mean(duration, na.rm=T),
					num_events=length(cel_sel),
					intensity_cumulative=sum(intensity_cumulative, na.rm=T),
					intensity_mean=mean(intensity_mean, na.rm=T),
					intensity_max=max(intensity_max, na.rm=T),
					intensity_var=mean(intensity_var, na.rm=T),
					rate_onset=mean(rate_onset, na.rm=T),
					rate_decline=mean(rate_decline, na.rm=T))
	out_df
	
}

rm(tda_dat)
rm(hw_dat_df)
gc()

outputName=paste('cov_', time_first, ".RData", sep="")
outputPath=file.path('~/tmp_cov_res', outputName)

assign(paste('cov_', time_first, sep=""), value=out, pos=1, inherits=T)
save(list=paste('cov_', time_first, sep=""), file=outputPath)

