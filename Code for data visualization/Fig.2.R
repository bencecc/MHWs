#### ---- Plots comparing mesoscale properties of MHWs networks ---- ####
#### ---- with null models ----------------------------------------- ####

# ----- Load libraries ------#

require(ggplot2)
require(dplyr)
require(ggsci)

# ----------- Load data with output from null models -----------#
setwd('~/MHWs/data')
load('tda_hist_null_stats.RData')
load('tda_hist_nulldep_stats.RData')
load('tda_rcp85_null_stats.RData')
load('tda_rcp85_nulldep_stats.RData')
load('tda_rcp26_null_stats.RData')
load('tda_rcp26_nulldep_stats.RData')
load('tda_glob_null_stats.RData')
load('tda_glob_nulldep_stats.RData')

# -------- Combine data ----------- #
# ---- Observed 
glob_null_stats <- tda_glob_null_stats
glob_null_stats_dep <- tda_glob_nulldep_stats

glob_obs_dat <- glob_null_stats[[1]][1:1000,]
glob_phase_dat <- glob_null_stats[[2]][1:1000,]
glob_dep_phase_dat <- glob_null_stats_dep[[2]][1:1000,]

glob_cond <- c(rep('random_phase', nrow(glob_phase_dat)),
		rep('constant_phase', nrow(glob_dep_phase_dat)),rep('observed', nrow(glob_obs_dat)))

glob_plot_null_dat <- data.frame(condition=glob_cond, rbind(glob_phase_dat,glob_dep_phase_dat,glob_obs_dat))

# ---- Historical
hist_null_stats <- tda_hist_null_stats
hist_null_stats_dep <- tda_hist_nulldep_stats

hist_obs_dat <- hist_null_stats[[1]][1:1000,]
hist_phase_dat <- hist_null_stats[[2]][1:1000,]
hist_dep_phase_dat <- hist_null_stats_dep[[2]][1:1000,]

hist_cond <- c(rep('random_phase', nrow(hist_phase_dat)),
		rep('constant_phase', nrow(hist_dep_phase_dat)), rep('observed', nrow(hist_obs_dat)))

hist_plot_null_dat <- data.frame(condition=hist_cond, rbind(hist_phase_dat,hist_dep_phase_dat,hist_obs_dat))

# ---- RCP2.6
rcp26_null_stats <- tda_rcp26_null_stats
rcp26_null_stats_dep <- tda_rcp26_nulldep_stats

rcp26_obs_dat <- rcp26_null_stats[[1]][1:1000,]
rcp26_phase_dat <- rcp26_null_stats[[2]][1:1000,]
rcp26_dep_phase_dat <- rcp26_null_stats_dep[[2]][1:1000,]

rcp26_cond <- c(rep('random_phase', nrow(rcp26_phase_dat)),
		rep('constant_phase', nrow(rcp26_dep_phase_dat)),rep('observed', nrow(rcp26_obs_dat)))

rcp26_plot_null_dat <- data.frame(condition=rcp26_cond, rbind(rcp26_phase_dat,rcp26_dep_phase_dat,rcp26_obs_dat))

# ---- RCP8.5
rcp85_null_stats <- tda_rcp85_null_stats
rcp85_null_stats_dep <- tda_rcp85_nulldep_stats

rcp85_obs_dat <- rcp85_null_stats[[1]][1:1000,]
rcp85_phase_dat <- rcp85_null_stats[[2]][1:1000,]
rcp85_dep_phase_dat <- rcp85_null_stats_dep[[2]][1:1000,]

rcp85_cond <- c(rep('random_phase', nrow(rcp85_phase_dat)),
		rep('constant_phase', nrow(rcp85_dep_phase_dat)), rep('observed', nrow(rcp85_obs_dat)))

rcp85_plot_null_dat <- data.frame(condition=rcp85_cond, rbind(rcp85_phase_dat,rcp85_dep_phase_dat,rcp85_obs_dat))

# ----- Merged dataset
null_models_dat <- data.frame(dataset=rep(c('observed','historical','rcp26','rcp85'), c(3000, 3000, 3000,3000)),
		rbind(glob_plot_null_dat,hist_plot_null_dat,rcp26_plot_null_dat,rcp85_plot_null_dat))

null_models_dat$dataset <- factor(null_models_dat$dataset)
null_models_dat$dataset <- factor(null_models_dat$dataset, levels(null_models_dat$dataset)[c(2,1,3,4)])
null_models_dat$condition <- factor(null_models_dat$condition)
null_models_dat$condition <- factor(null_models_dat$condition, levels(null_models_dat$condition)[c(3,1,2)])

int.x=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5,6.5,7,7.5)
reps <- c(1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000)
null_models_dat <- data.frame(int_x=rep(int.x, reps), null_models_dat)

code <- rep(c('obs_rand','obs_const','obs_obs','hist_rand','hist_const','hist_obs',
				'rcp26_rand','rcp26_const','rcp26_obs','rcp85_rand','rcp85_const','rcp85_obs'), reps)

null_models_dat <- data.frame(cbind(code, null_models_dat))
null_models_dat$code <- factor(null_models_dat$code)
null_models_dat$code <- factor(null_models_dat$code, levels(null_models_dat$code)[c(6,4,5,3,1,2,9,7,8,12,10,11)])

obs_means <- null_models_dat %>% filter(condition=='observed') %>% group_by(int_x, code, dataset, condition) %>%
		summarise_all(mean)

int.x=c(0.5,1,1.5,2.5,3,3.5,4.5,5,5.5,6.5,7,7.5)

p_mod_null <- ggplot(data=null_models_dat %>% filter(condition!='observed'), aes(x=int_x, y=modularity, group=code))+
		geom_violin(trim=F, show.legend=F, scale='width', aes(fill=code))+
		geom_boxplot(width=0.2, show.legend=F, aes(fill=code))+
		geom_point(data=obs_means, aes(x=int_x, y=modularity, fill=dataset), shape=23, size=5, show.legend=F)+
		scale_x_continuous(name="", breaks=int.x, labels=rep(c('RP','CP','OB'), 4))+
		theme_bw()+
		theme(axis.text.x=element_text(angle=0, size=11),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank())

windows(width=5, height=4)
p_mod_null

#ggsave("~/MHWs/data/modularity_null.pdf", useDingbats=FALSE, width = 5, height = 4)

p_deg_null <- ggplot(data=null_models_dat %>% filter(condition!='observed'), aes(x=int_x, y=mean_deg, group=code))+
		geom_violin(trim=F, show.legend=F, scale='width', aes(fill=code))+
		geom_boxplot(width=0.2, show.legend=F, aes(fill=code))+
		geom_point(data=obs_means, aes(x=int_x, y=mean_deg, fill=code), shape=23, size=5, show.legend=F)+
		scale_x_continuous(name="", breaks=int.x, labels=rep(c('RP','CP','OB'), 4))+
		theme_bw()+
		theme(axis.text.x=element_text(angle=0, size=11),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank())

windows(width=6, height=5)
p_deg_null

#ggsave("~/MHWs/data/degree_null.pdf", useDingbats=FALSE, width = 6, height = 5)


# --------- Signficance tests ----------- #
# These are two-tailed tests: Exc= n>=(<=)OBS; 2*Exc/n;
# in addition, Bonferroni correction was employed to account for
# multiple testing (prob = 2*outcome two-tailed test) 
# modularity
pvals_mod_obs_rand <- null_models_dat %>% filter(condition!='observed'&dataset=='observed') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase','modularity']>=
										as.numeric(obs_means[obs_means$code=='obs_obs','modularity'])))/1000))
pvals_mod_obs_const <- null_models_dat %>% filter(condition!='observed'&dataset=='observed') %>%
		do(data.frame(p_value=length(which(.[.$condition=='constant_phase','modularity']>=
										as.numeric(obs_means[obs_means$code=='obs_obs','modularity'])))/1000))
pvals_mod_hist_rand <- null_models_dat %>% filter(condition!='observed'&dataset=='historical') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase','modularity']<=
										as.numeric(obs_means[obs_means$code=='hist_obs','modularity'])))/1000))
pvals_mod_hist_const <- null_models_dat %>% filter(condition!='observed'&dataset=='historical') %>%
		do(data.frame(p_value=length(which(.[.$condition=='constant_phase','modularity']<=
										as.numeric(obs_means[obs_means$code=='hist_obs','modularity'])))/1000))
pvals_mod_rcp26_rand <- null_models_dat %>% filter(condition!='observed'&dataset=='rcp26') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase','modularity']>=
										as.numeric(obs_means[obs_means$code=='rcp26_obs','modularity'])))/1000))
pvals_mod_rcp26_const <- null_models_dat %>% filter(condition!='observed'&dataset=='rcp26') %>%
		do(data.frame(p_value=length(which(.[.$condition=='constant_phase','modularity']>=
										as.numeric(obs_means[obs_means$code=='rcp26_obs','modularity'])))/1000))
pvals_mod_rcp85_rand <- null_models_dat %>% filter(condition!='observed'&dataset=='rcp85') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase','modularity']<=
												as.numeric(obs_means[obs_means$code=='rcp85_obs','modularity'])))/1000))
pvals_mod_rcp85_const <- null_models_dat %>% filter(condition!='observed'&dataset=='rcp85') %>%
		do(data.frame(p_value=length(which(.[.$condition=='constant_phase','modularity']<=
												as.numeric(obs_means[obs_means$code=='rcp85_obs','modularity'])))/1000))

# degree
pvals_deg_obs_rand <- null_models_dat %>% filter(condition!='observed'&dataset=='observed') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase','mean_deg']<=
												as.numeric(obs_means[obs_means$code=='obs_obs','mean_deg'])))/1000))
pvals_deg_obs_const <- null_models_dat %>% filter(condition!='observed'&dataset=='observed') %>%
		do(data.frame(p_value=length(which(.[.$condition=='constant_phase','mean_deg']<=
												as.numeric(obs_means[obs_means$code=='obs_obs','mean_deg'])))/1000))
pvals_deg_hist_rand <- null_models_dat %>% filter(condition!='observed'&dataset=='historical') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase','mean_deg']>=
												as.numeric(obs_means[obs_means$code=='hist_obs','mean_deg'])))/1000))
pvals_deg_hist_const <- null_models_dat %>% filter(condition!='observed'&dataset=='historical') %>%
		do(data.frame(p_value=length(which(.[.$condition=='constant_phase','mean_deg']>=
												as.numeric(obs_means[obs_means$code=='hist_obs','mean_deg'])))/1000))
pvals_deg_rcp26_rand <- null_models_dat %>% filter(condition!='observed'&dataset=='rcp26') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase','mean_deg']<=
												as.numeric(obs_means[obs_means$code=='rcp26_obs','mean_deg'])))/1000))
pvals_deg_rcp26_const <- null_models_dat %>% filter(condition!='observed'&dataset=='rcp26') %>%
		do(data.frame(p_value=length(which(.[.$condition=='constant_phase','mean_deg']<=
												as.numeric(obs_means[obs_means$code=='rcp26_obs','mean_deg'])))/1000))
pvals_deg_rcp85_rand <- null_models_dat %>% filter(condition!='observed'&dataset=='rcp85') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase','mean_deg']>=
												as.numeric(obs_means[obs_means$code=='rcp85_obs','mean_deg'])))/1000))
pvals_deg_rcp85_const <- null_models_dat %>% filter(condition!='observed'&dataset=='rcp85') %>%
		do(data.frame(p_value=length(which(.[.$condition=='constant_phase','mean_deg']>=
												as.numeric(obs_means[obs_means$code=='rcp85_obs','mean_deg'])))/1000))

# --------------------------------------- #


