#### ----------- Reproduce Supplementary Figs. 5-8. Plots of network graphs ----- ####
#### ----------- under different gain and resolution parameters ----------------- ####

# ---- Require libraries
require(abind)
require(dplyr)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)
require(ggsci)
library(wesanderson)
require(viridis)
require(TDAmapper)
require(igraph)
require(reticulate)

setwd('~/MHWs/data')

# ---- Load functions
source('~/MHWs/code/plot_net.R') #function to generate network plots from TDA objects
source('~/MHWs/code/color_nodes_time.R') # function to color nodes in the network
source('~/MHWs/code/net_stats.R') # function to calculate network statistics
source('~/MHWs/code/temporal_degree.R') # function to compute Normalized degree from a temporal connectivity matrix (TCM)

# ---- Observed SST
# covariate data to color nodes
load('~/MHWs/data/tda_glob_covar.RData')
# load objects returned by the TDA mapper algorithm for several combinations of resolution and gain 
load('~/MHWs/data/glob_net_pertur_dat.RData')

windows(height=14, width=14)
ggarrange(plotlist=glob_net_pertur_dat,
		ncol=4, nrow=4,
		common.legend=T,
		legend=TRUE
)

#ggsave("~/MHWs/data/glob_perturb_net_tmp.pdf", useDingbats=FALSE, width = 14, height = 14)

# ---- Statistics and data on extreme networks
load('~/MHWs/data/tda_glob_null_stats.RData')
load('~/MHWs/data/glob_years_length.RData')
load('~/MHWs/data/glob_res_mat.RData')
load('~/MHWs/data/glob_res_mat.RData') 
load('~/MHWs/data/tda_glob_12_25.RData')
load('~/MHWs/data/tda_glob_30_55.RData')

norm_deg_glob <- glob_res_mat[[2]]
glob1 <- tda_glob_12_25
glob2 <- tda_glob_30_55

tmp <- temporal_degree(glob1) # function temporal_degree uses requires foreach and doMC
node_degree_glob1 <- (tmp-min(tmp))/(max(tmp)-min(tmp))
tmp1 <- temporal_degree(glob2)
node_degree_glob2 <- (tmp1-min(tmp1))/(max(tmp1)-min(tmp1))
comp_glob_deg <- data.frame(time=rep(1:length(norm_deg_glob),3), net_type=rep(c("12_25","24_45","30_55"), each=length(norm_deg_glob)),
		deg=c(node_degree_glob1[1:13513],norm_deg_glob, node_degree_glob2[1:13513]))

# null data
glob_df <- as.data.frame(tda_glob_null_stats[[4]])
glob_df_long <- stack(glob_df)
glob_deg_df <- data.frame(time=rep(1:13514, nrow(glob_df)), model=rep(1:nrow(glob_df), each=13514), degree=glob_df_long$value) 

deg_tmp <- apply(glob_df, 2, function(x) {
			mean=mean(x)
			sd=sd(x)
			c(mean, sd)
		})

deg_tmp <- t(deg_tmp)
deg_null <- data.frame(time=1:nrow(deg_tmp), mean_deg=deg_tmp[,1],
		sd=deg_tmp[,2])

plot_years <- c(1982,1985,1990,1995,2000,2005,2010,2015,2018)
int.time <- sapply(plot_years, function(x) {
			sel <- which(glob_years_length$year%in%x)
			sum(years_length_obs[1:sel,'no_days'])	
		})

degree_glob_p <- ggplot(deg_null, aes(x=time, y=mean_deg))+
		geom_line(data=comp_glob_deg, aes(x=time, y=deg, color=factor(net_type)), alpha=0.2, linetype=4)+
		stat_smooth(data=comp_glob_deg, aes(x=time, y=deg, color=factor(net_type)), se=T)+
		geom_ribbon(aes(ymin=ifelse(mean_deg-2*sd>=0,(mean_deg-2*sd),0),
					ymax=mean_deg+2*sd), fill="#ACD3A2", alpha=0.3)+
		geom_line(colour='white')+
		scale_color_manual(name = "Resolution \nand overlap", values = c("darkorange2", 
						"firebrick", 
						"dodgerblue3"))+
		scale_x_continuous(name="", breaks=c(int.time[2:8]), expand=c(0.05,0.05), 
				labels=c('1985','1990','1995','2000','2005','2010','2015'))+
		scale_y_continuous(name="Normalized degree of TCM", limits=c(-0.05,1.01), expand=c(0,0))+
		theme_bw(base_size = 12)+
		theme(panel.grid.major=element_line(colour=NA),
				legend.title=element_text(size=10),
				axis.text.x=element_text(size=11),
				axis.text.y=element_text(size=11),
				axis.title.y=element_text(size=11),
				panel.grid.minor=element_line(colour=NA))

windows(width=6, height=3)
degree_glob_p

#ggsave("~/MHWs/data/glob_perturb_stats.pdf", useDingbats=FALSE, width = 6, height = 3)


# ---- Historical
# load covariate data to color nodes
load('~MHWs/data/tda_hist_covar.RData')
# load objects returned by the TDA mapper algorithm for many combinations of resolution and gain 
load('~/MHWs/data/hist_net_pertur_dat.RData')
	
windows(height=14, width=14)
	ggarrange(plotlist=hist_net_pertur_dat,
			ncol=4, nrow=4,
			common.legend=T,
			legend=T
			)

#ggsave("~/MHWs/data/hist_perturb_net_tmp.pdf", useDingbats=FALSE, width = 14, height = 14)

# ---- Statistics on extreme networks
load('~/MHWs/data/tda_hist_null_stats.RData')
load('~/MHWs/data/hist_years_length.RData')
load('~/MHWs/data/hist_res_mat.RData') #load('glob_res_hist_mat.RData')
load('~/MHWs/data/tda_hist_10_25.RData')
load('~/MHWs/data/tda_hist_28_55.RData')

norm_deg_hist <- hist_res_mat[[2]]
hist1 <- tda_hist_10_25
hist2 <- tda_hist_28_55

tmp <- temporal_degree(hist1)
node_degree_hist1 <- (tmp-min(tmp))/(max(tmp)-min(tmp))
tmp1 <- temporal_degree(hist2)
node_degree_hist2 <- (tmp1-min(tmp1))/(max(tmp1)-min(tmp1))
comp_hist_deg <- data.frame(time=rep(1:52959,3), net_type=rep(c("10_25","30_55","28_55"), each=52959),
		deg=c(node_degree_hist1[1:52959],norm_deg_hist[1:52959], node_degree_hist2))

# null data
hist_df <- as.data.frame(tda_hist_null_stats[[4]])
hist_df_long <- stack(hist_df)
hist_deg_df <- data.frame(time=rep(1:52960, nrow(hist_df)), model=rep(1:nrow(hist_df), each=52960), degree=hist_df_long$value) 

deg_tmp <- apply(hist_df, 2, function(x) {
			mean=mean(x)
			sd=sd(x)
			c(mean, sd)
		})

deg_tmp <- t(deg_tmp)
deg_null <- data.frame(time=1:nrow(deg_tmp), mean_deg=deg_tmp[,1],
		sd=deg_tmp[,2])

plot_years <- c(1880,1900,1920,1940,1960,1980,2000)
int.time <- sapply(plot_years, function(x) {
			sel <- which(hist_years_length$year%in%x)
			sum(hist_years_length[1:sel,'no_days'])	
		})

degree_hist_p <- ggplot(deg_null, aes(x=time, y=mean_deg))+
		geom_line(data=comp_hist_deg, aes(x=time, y=deg, color=factor(net_type)), alpha=0.2, linetype=4)+
		stat_smooth(data=comp_hist_deg, aes(x=time, y=deg, color=factor(net_type)), se=T)+
		geom_ribbon(aes(ymin=ifelse(mean_deg-2*sd>=0,(mean_deg-2*sd),0),
						ymax=mean_deg+2*sd), fill="#ACD3A2", alpha=0.5)+
		geom_line(colour='white')+
		scale_color_manual(name = "Resolution \nand overlap", values = c("darkorange2", 
						"firebrick", 
						"dodgerblue3"))+
		scale_x_continuous(name="", breaks=c(1,int.time), expand=c(0.05,0.05), 
				labels=c('1861','1880','1900','1920','1940','1960','1980','2000'))+
		scale_y_continuous(name="Normalized degree of TCM", limits=c(-0.05,1.01), expand=c(0,0))+
		theme_bw(base_size = 12)+
		theme(panel.grid.major=element_line(colour=NA),
				legend.title=element_text(size=10),
				axis.text.x=element_text(size=11),
				axis.text.y=element_text(size=11),
				axis.title.y=element_text(size=11),
				panel.grid.minor=element_line(colour=NA))


windows(width=6, height=3)
degree_hist_p

#ggsave("~/MHWs/data/hist_perturb_stats.pdf", useDingbats=FALSE, width = 6, height = 3)

# ---- RCP 8.5
# covariate data to color nodes
load('~/MHWs/data/tda_rcp85_covar.RData')
load('~/MHWs/data/rcp85_net_perturb_dat.RData')

windows(height=14, width=14)
ggarrange(plotlist=rcp85_net_perturb_dat,
		ncol=4, nrow=4,
		common.legend=T
#		legend=FALSE
)

#ggsave("~/MHWs/data/rcp85_perturb_net_legend.pdf", useDingbats=FALSE, width = 14, height = 14)

# ---- Statistics on extreme networks
load('~/MHWs/data/tda_rcp85_null_stats.RData')
load('~/MHWs/data/rcp85_years_length.RData')
load('~/MHWs/data/rcp85_res_mat.RData')
load('~/MHWs/data/tda_rcp85_12_25.RData')
load('~/MHWs/data/tda_rcp85_30_55.RData')

norm_deg_rcp85 <- rcp85_res_mat[[2]]
rcp851 <- tda_rcp85_12_25
rcp852 <- tda_rcp85_30_55

stat_rcp851 <- net_stats(rcp851)
stat_rcp852 <- net_stats(rcp852)

tmp <- temporal_degree(rcp851)
node_degree_rcp851 <- (tmp-min(tmp))/(max(tmp)-min(tmp))
tmp1 <- temporal_degree(rcp852)
node_degree_rcp852 <- (tmp1-min(tmp1))/(max(tmp1)-min(tmp1))
comp_rcp85_deg <- data.frame(time=rep(1:34674,3), net_type=rep(c("12_25","24_45","30_55"), each=34674),
		deg=c(node_degree_rcp851[1:34674],norm_deg_rcp85[1:34674], node_degree_rcp852))

# null data
rcp85_df <- as.data.frame(tda_rcp85_null_stats[[4]])
rcp85_df_long <- stack(rcp85_df)
rcp85_deg_df <- data.frame(time=rep(1:34675, nrow(rcp85_df)), model=rep(1:nrow(rcp85_df), each=34675), degree=rcp85_df_long$value) 

deg_tmp <- apply(rcp85_df, 2, function(x) {
			mean=mean(x)
			sd=sd(x)
			c(mean, sd)
		})

deg_tmp <- t(deg_tmp)
deg_null <- data.frame(time=1:nrow(deg_tmp), mean_deg=deg_tmp[,1],
		sd=deg_tmp[,2])

plot_years <- c(2010,2020,2030,2040,2050,2060,2070,2080,2090,2100)
int.time <- sapply(plot_years, function(x) {
			sel <- which(rcp85_years_length$year%in%x)
			sum(rcp85_years_length[1:sel,'no_days'])	
		})

degree_rcp85_p <- ggplot(deg_null, aes(x=time, y=mean_deg))+
		geom_line(data=comp_rcp85_deg, aes(x=time, y=deg, color=factor(net_type)), alpha=0.2, linetype=4)+
		stat_smooth(data=comp_rcp85_deg, aes(x=time, y=deg, color=factor(net_type)), se=T)+
		geom_ribbon(aes(ymin=ifelse(mean_deg-2*sd>=0,(mean_deg-2*sd),0),
						ymax=mean_deg+2*sd), fill="#ACD3A2", alpha=0.5)+
		geom_line(colour='white')+
		scale_color_manual(name = "Resolution \nand overlap", values = c("darkorange2", 
						"firebrick", 
						"dodgerblue3"))+
		scale_x_continuous(name="", breaks=int.time, expand=c(0,0), 
				labels=c('2010','2020','2030','2040','2050','2060','2070',
						'2080','2090','2100'))+
		scale_y_continuous(name="Normalized degree of TCM", limits=c(-0.05,1.01), expand=c(0,0))+
		theme_bw(base_size = 12)+
		theme(panel.grid.major=element_line(colour=NA),
				legend.title=element_text(size=10),
				axis.text.x=element_text(size=11),
				axis.text.y=element_text(size=11),
				axis.title.y=element_text(size=11),
				panel.grid.minor=element_line(colour=NA))


windows(width=6, height=3)
degree_rcp85_p

#ggsave("~/MHWs/data/rcp85_perturb_stats.pdf", useDingbats=FALSE, width = 6, height = 3)

# ---- RCP 2.6
# covariate data to color nodes
load('~/MHWs/data/tda_rcp26_covar.RData')
load('rcp26_net_pertur_dat.RData')

windows(height=14, width=14)
ggarrange(plotlist=rcp26_net_pertur_dat,
		ncol=4, nrow=4,
		#common.legend=T,
		legend=FALSE
)

#ggsave("~/MHWs/data/rcp26_perturb_net_tmp.pdf", useDingbats=FALSE, width = 14, height = 14)

# ---- Statistics on extreme networks
load('~/MHWs/data/tda_rcp26_null_stats.RData')
load('~/MHWs/data/rcp26_years_length.RData')
load('~/MHWs/data/rcp26_res_mat.RData')
load('~/MHWs/data/tda_rcp26_6_15.RData')
load('~/MHWs/data/tda_rcp26_24_45.RData')

norm_deg_rcp26 <- rcp26_res_mat[[2]]
rcp261 <- tda_rcp26_6_15
rcp262 <- tda_rcp26_24_45

tmp <- temporal_degree(rcp261)
node_degree_rcp261 <- (tmp-min(tmp))/(max(tmp)-min(tmp))
tmp1 <- temporal_degree(rcp262)
node_degree_rcp262 <- (tmp1-min(tmp1))/(max(tmp1)-min(tmp1))
comp_rcp26_deg <- data.frame(time=rep(1:34673,3), net_type=rep(c("6_15","12_25","24_45"), each=34673),
		deg=c(node_degree_rcp261[1:34673],norm_deg_rcp26[1:34673], node_degree_rcp262[1:34673]))
comp_rcp26_deg$net_type <- factor(comp_rcp26_deg$net_type, levels=c("6_15","12_25","24_45"))

# null data
rcp26_df <- as.data.frame(tda_rcp26_null_stats[[4]])
rcp26_df_long <- stack(rcp26_df)
rcp26_deg_df <- data.frame(time=rep(1:34675, nrow(rcp26_df)), model=rep(1:nrow(rcp26_df), each=34675),
		degree=rcp26_df_long$value) 

deg_tmp <- apply(rcp26_df, 2, function(x) {
			mean=mean(x)
			sd=sd(x)
			c(mean, sd)
		})

deg_tmp <- t(deg_tmp)
deg_null <- data.frame(time=1:nrow(deg_tmp), mean_deg=deg_tmp[,1],
		sd=deg_tmp[,2])

plot_years <- c(2010,2020,2030,2040,2050,2060,2070,2080,2090,2100)
int.time <- sapply(plot_years, function(x) {
			sel <- which(rcp26_years_length$year%in%x)
			sum(rcp26_years_length[1:sel,'no_days'])	
		})

degree_rcp26_p <- ggplot(deg_null, aes(x=time, y=mean_deg))+
		geom_line(data=comp_rcp26_deg, aes(x=time, y=deg, color=factor(net_type)), alpha=0.2, linetype=4)+
		stat_smooth(data=comp_rcp26_deg, aes(x=time, y=deg, color=factor(net_type)), se=T)+
		geom_ribbon(aes(ymin=ifelse(mean_deg-2*sd>=0,(mean_deg-2*sd),0),
						ymax=mean_deg+2*sd), fill="#ACD3A2", alpha=0.5)+
		geom_line(colour='white')+
		scale_color_manual(name = "Resolution \nand overlap", values = c("darkorange2", 
						"firebrick", 
						"dodgerblue3"))+
		scale_x_continuous(name="", breaks=int.time, expand=c(0,0), 
				labels=c('2010','2020','2030','2040','2050','2060','2070',
						'2080','2090','2100'))+
		scale_y_continuous(name="Normalized degree of TCM", limits=c(-0.05,1.01), expand=c(0,0))+
		theme_bw(base_size = 12)+
		theme(panel.grid.major=element_line(colour=NA),
				legend.title=element_text(size=10),
				axis.text.x=element_text(size=11),
				axis.text.y=element_text(size=11),
				axis.title.y=element_text(size=11),
				panel.grid.minor=element_line(colour=NA))


windows(width=6, height=3)
degree_rcp26_p

#ggsave("~/MHWs/data/rcp26_perturb_stats.pdf", useDingbats=FALSE, width = 6, height = 3)

#### ---- NULL MODEL PLOTS ------------------------------------------ ####

# ----- Observed SST 
load('~/MHWs/data/tda_glob_null_stats.RData')
load('~/MHWs/data/glob_optim1225_null_stats.RData')
load('~/MHWs/data/glob_optim3055_null_stats.RData')
		
glob_null_stats <- tda_glob_null_stats
glob_obs_dat <- glob_null_stats[[1]][1:1000,]
glob_phase_dat <- glob_null_stats[[2]][1:1000,]

glob_obs_optim1225 <- glob_optim1225_null_stats[[1]][1:1000,]
glob_phase_optim1225 <- glob_optim1225_null_stats[[2]][1:1000,]
 
glob_obs_optim3055 <- glob_optim3055_null_stats[[1]][1:1000,]
glob_phase_optim3055 <- glob_optim3055_null_stats[[2]][1:1000,]


glob_cond <- c(rep('random_phase_12_25', nrow(glob_phase_optim1225)),rep('observed_12_25', nrow(glob_phase_optim1225)),
		rep('random_phase_24_45', nrow(glob_phase_dat)), rep('observed_24_45', nrow(glob_obs_dat)),
		rep('random_phase_30_55', nrow(glob_phase_optim3055)), rep('observed_30_55', nrow(glob_phase_optim3055)))

int.x=c(0,1,3,4,6,7)
glob_plot_null_dat <- data.frame(int_x=rep(c(0,1,3,4,6,7), each=1000), condition=glob_cond, rbind(glob_phase_optim1225,glob_obs_optim1225,
				glob_phase_dat,glob_obs_dat, glob_phase_optim3055, glob_obs_optim3055))

glob_plot_null_dat$condition <- factor(glob_plot_null_dat$condition)
glob_plot_null_dat$condition <- factor(glob_plot_null_dat$condition, levels(glob_plot_null_dat$condition)[c(4,1,5,2,6,3)])

obs_vals <- glob_plot_null_dat %>% filter(condition==c('observed_12_25','observed_24_45','observed_30_55')) %>%
		group_by(condition) %>% summarise_all(mean)

glob_mod_null_perturbed <- ggplot(data=glob_plot_null_dat %>%
						filter(condition==c('random_phase_12_25','random_phase_24_45', 'random_phase_30_55')),
				aes(x=int_x, y=modularity, group=condition))+
		geom_violin(trim=T, show.legend=F, scale='width', aes(fill=condition), color='black', size=0.3)+
		geom_boxplot(width=0.2, show.legend=F, aes(fill=condition))+
		geom_point(data=obs_vals, aes(x=int_x, y=modularity, fill=condition), color='black', shape=23, size=5, show.legend=F)+
		scale_fill_manual(name = "", values = c("darkorange2", "firebrick", "dodgerblue3",
						"darkorange2","firebrick", "dodgerblue3"))+
		scale_color_manual(name = "", values = c("darkorange2", "firebrick", "dodgerblue3",
						"darkorange2","firebrick", "dodgerblue3"))+
		scale_x_continuous(name="", breaks=int.x,
				labels=c('RP\n12_25','OB\n12_25','RP\n24_45','OB\n24_45','RP\n30_55','OB\n30_55'),
				expand=c(0.05,0.05))+
		scale_y_continuous(name="Modularity", limits=c(0.87,1))+
		theme_bw()+
		theme(axis.text.x=element_text(angle=0, size=8),
				axis.text.y=element_text(size=10),
				axis.title.y=element_text(size=10),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank())

windows(width=4, height=5)
glob_mod_null_perturbed

#ggsave("~/MHWs/data/glob_perturb_modularity.pdf", useDingbats=FALSE, width = 4, height = 5)

pvals_mod_12_25 <- glob_plot_null_dat %>% filter(condition=='random_phase_12_25') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_12_25','modularity']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_12_25','modularity'])))/1000))
pvals_mod_24_45 <- glob_plot_null_dat %>% filter(condition=='random_phase_24_45') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_24_45','modularity']>=
												as.numeric(obs_vals[obs_vals$condition=='observed_24_45','modularity'])))/1000))
pvals_mod_30_55 <- glob_plot_null_dat %>% filter(condition=='random_phase_30_55') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_30_55','modularity']>=
												as.numeric(obs_vals[obs_vals$condition=='observed_30_55','modularity'])))/1000))

glob_deg_null_perturbed <- ggplot(data=glob_plot_null_dat %>%
						filter(condition==c('random_phase_12_25','random_phase_24_45', 'random_phase_30_55')),
				aes(x=int_x, y=mean_deg, group=condition))+
		geom_violin(trim=T, show.legend=F, scale='width', aes(fill=condition), color='black', size=0.3)+
		geom_boxplot(width=0.2, show.legend=F, aes(fill=condition))+
		geom_point(data=obs_vals, aes(x=int_x, y=mean_deg, fill=condition), color='black', shape=23, size=5, show.legend=F)+
		scale_fill_manual(name = "", values = c("darkorange2", "firebrick", "dodgerblue3",
						"darkorange2","firebrick", "dodgerblue3"))+
		scale_color_manual(name = "", values = c("darkorange2", "firebrick", "dodgerblue3",
						"darkorange2","firebrick", "dodgerblue3"))+
		scale_x_continuous(name="", breaks=int.x,
				labels=c('RP\n12_25','OB\n12_25','RP\n24_45','OB\n24_45','RP\n30_55','OB\n30_55'),
				expand=c(0.05,0.05))+
		scale_y_continuous(name="Node degree", limits=c(1.8,3.2))+
		theme_bw()+
		theme(axis.text.x=element_text(angle=0, size=8),
				axis.text.y=element_text(size=10),
				axis.title.y=element_text(size=10),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank())

windows(width=4, height=5)
glob_deg_null_perturbed		

#ggsave("~/MHWs/data/glob_perturb_degree.pdf", useDingbats=FALSE, width = 4, height = 5)

pvals_deg_12_25 <- glob_plot_null_dat %>% filter(condition=='random_phase_12_25') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_12_25','mean_deg']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_12_25','mean_deg'])))/1000))
pvals_deg_24_45 <- glob_plot_null_dat %>% filter(condition=='random_phase_24_45') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_24_45','mean_deg']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_24_45','mean_deg'])))/1000))
pvals_deg_30_55 <- glob_plot_null_dat %>% filter(condition=='random_phase_30_55') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_30_55','mean_deg']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_30_55','mean_deg'])))/1000))

		
# ----- Historical 
load('~/MHWs/data/tda_hist_null_stats.RData')
load('~/MHWs/data/hist_optim1025_null_stats.RData')
load('~/MHWs/data/hist_optim2855_null_stats.RData')

hist_null_stats <- tda_hist_null_stats
hist_obs_dat <- hist_null_stats[[1]][1:1000,]
hist_phase_dat <- hist_null_stats[[2]][1:1000,]

hist_obs_optim1025 <- hist_optim1025_null_stats[[1]][1:1000,]
hist_phase_optim1025 <- hist_optim1025_null_stats[[2]][1:1000,]

hist_obs_optim2855 <- hist_optim2855_null_stats[[1]][1:1000,]
hist_phase_optim2855 <- hist_optim2855_null_stats[[2]][1:1000,]

hist_cond <- c(rep('random_phase_10_25', nrow(hist_phase_optim1025)),rep('observed_10_25', nrow(hist_phase_optim1025)),
		rep('random_phase_30_55', nrow(hist_phase_dat)), rep('observed_30_55', nrow(hist_obs_dat)),
		rep('random_phase_28_55', nrow(hist_phase_optim2855)), rep('observed_28_55', nrow(hist_phase_optim2855)))

int.x=c(0,1,3,4,6,7)
hist_plot_null_dat <- data.frame(int_x=c(rep(c(0,1), each=1000), rep(c(3,4), each=1000), rep(c(6,7), each=1000)),
		condition=hist_cond, rbind(hist_phase_optim1025,hist_obs_optim1025,
				hist_phase_dat,hist_obs_dat, hist_phase_optim2855, hist_obs_optim2855))

hist_plot_null_dat$condition <- factor(hist_plot_null_dat$condition)
hist_plot_null_dat$condition <- factor(hist_plot_null_dat$condition, levels(hist_plot_null_dat$condition)[c(4,1,5,2,6,3)])

obs_vals <- hist_plot_null_dat %>% filter(condition==c('observed_10_25','observed_30_55','observed_28_55')) %>%
		group_by(condition) %>% summarise_all(mean)

hist_mod_null_perturbed <- ggplot(data=hist_plot_null_dat %>%
						filter(condition==c('random_phase_10_25','random_phase_30_55', 'random_phase_28_55')),
				aes(x=int_x, y=modularity))+
		geom_violin(trim=T, show.legend=F, scale='width', aes(fill=condition), color='black', size=0.3)+
		geom_boxplot(width=0.2, show.legend=F, aes(fill=condition))+
		geom_point(data=obs_vals, aes(x=int_x, y=modularity, fill=condition), color='black', shape=23, size=5, show.legend=F)+
		scale_fill_manual(name = "", values = c("darkorange2", "dodgerblue3","firebrick", 
						"darkorange2","dodgerblue3","firebrick"))+
		scale_color_manual(name = "", values = c("darkorange2", "dodgerblue3","firebrick", 
						"darkorange2", "dodgerblue3", "firebrick"))+
		scale_x_continuous(name="", breaks=int.x,
				labels=c('RP\n10_25','OB\n10_25','RP\n30_55','OB\n30_55','RP\n28_55','OB\n28_55'),
				expand=c(0.05,0.05))+
		scale_y_continuous(name="Modularity", limits=c(0.76,0.98))+
		theme_bw()+
		theme(axis.text.x=element_text(angle=0, size=8),
				axis.text.y=element_text(size=10),
				axis.title.y=element_text(size=10),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank())

windows(width=4, height=5)
hist_mod_null_perturbed

#ggsave("~/MHWs/data/hist_perturb_modularity.pdf", useDingbats=FALSE, width = 4, height = 5)

pvals_mod_10_25 <- hist_plot_null_dat %>% filter(condition=='random_phase_10_25') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_10_25','modularity']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_10_25','modularity'])))/1000))
pvals_mod_30_55 <- hist_plot_null_dat %>% filter(condition=='random_phase_30_55') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_30_55','modularity']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_30_55','modularity'])))/1000))
pvals_mod_28_55 <- hist_plot_null_dat %>% filter(condition=='random_phase_28_55') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_28_55','modularity']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_28_55','modularity'])))/1000))

hist_deg_null_perturbed <- ggplot(data=hist_plot_null_dat %>%
						filter(condition==c('random_phase_10_25','random_phase_30_55', 'random_phase_28_55')),
				aes(x=int_x, y=mean_deg, group=condition))+
		geom_violin(trim=T, show.legend=F, scale='width', aes(fill=condition), color='black')+
		geom_boxplot(width=0.2, show.legend=F, aes(fill=condition))+
		geom_point(data=obs_vals, aes(x=int_x, y=mean_deg, fill=condition), color='black', shape=23, size=5, show.legend=F)+
		scale_fill_manual(name = "", values = c("darkorange2", "dodgerblue3","firebrick", 
						"darkorange2","dodgerblue3","firebrick"))+
		scale_color_manual(name = "", values = c("darkorange2", "dodgerblue3","firebrick", 
						"darkorange2", "dodgerblue3", "firebrick"))+
		scale_x_continuous(name="", breaks=int.x,
				labels=c('RP\n10_25','OB\n10_25','RP\n30_55','OB\n30_55','RP\n28_55','OB\n28_55'),
				expand=c(0.05,0.05))+
		scale_y_continuous(name="Node degree")+
		theme_bw()+
		theme(axis.text.x=element_text(angle=0, size=8),
				axis.text.y=element_text(size=10),
				axis.title.y=element_text(size=10),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank())

windows(width=4, height=5)
hist_deg_null_perturbed

#ggsave("~/MHWs/data/hist_perturb_degree.pdf", useDingbats=FALSE, width = 4, height = 5)

pvals_deg_10_25 <- hist_plot_null_dat %>% filter(condition=='random_phase_10_25') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_10_25','mean_deg']>=
												as.numeric(obs_vals[obs_vals$condition=='observed_10_25','mean_deg'])))/1000))
pvals_deg_30_55 <- hist_plot_null_dat %>% filter(condition=='random_phase_30_55') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_30_55','mean_deg']>=
												as.numeric(obs_vals[obs_vals$condition=='observed_30_55','mean_deg'])))/1000))
pvals_deg_28_55 <- hist_plot_null_dat %>% filter(condition=='random_phase_28_55') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_28_55','mean_deg']>=
												as.numeric(obs_vals[obs_vals$condition=='observed_28_55','mean_deg'])))/1000))

# ----- RCP 2.6 
load('~/MHWs/data/tda_rcp26_null_stats.RData')
load('~/MHWs/data/rcp26_optim615_null_stats.RData')
load('~/MHWs/data/rcp26_optim2445_null_stats.RData')

rcp26_null_stats <- tda_rcp26_null_stats
rcp26_obs_dat <- rcp26_null_stats[[1]][1:1000,]
rcp26_phase_dat <- rcp26_null_stats[[2]][1:1000,]

rcp26_obs_optim615 <- rcp26_optim615_null_stats[[1]][1:1000,]
rcp26_phase_optim615 <- rcp26_optim615_null_stats[[2]][1:1000,]

rcp26_obs_optim2445 <- rcp26_optim2445_null_stats[[1]][1:1000,]
rcp26_phase_optim2445 <- rcp26_optim2445_null_stats[[2]][1:1000,]

rcp26_cond <- c(rep('random_phase_6_15', nrow(rcp26_phase_optim615)),rep('observed_6_15', nrow(rcp26_phase_optim615)),
		rep('random_phase_12_25', nrow(rcp26_phase_dat)), rep('observed_12_25', nrow(rcp26_obs_dat)),
		rep('random_phase_24_45', nrow(rcp26_phase_optim2445)), rep('observed_24_45', nrow(rcp26_phase_optim2445)))

int.x=c(0,1,3,4,6,7)
rcp26_plot_null_dat <- data.frame(int_x=c(rep(c(0,1), each=1000), rep(c(3,4), each=1000), rep(c(6,7), each=1000)),
		condition=rcp26_cond, rbind(rcp26_phase_optim615,rcp26_obs_optim615,
				rcp26_phase_dat,rcp26_obs_dat, rcp26_phase_optim2445, rcp26_obs_optim2445))

rcp26_plot_null_dat$condition <- factor(rcp26_plot_null_dat$condition)
rcp26_plot_null_dat$condition <- factor(rcp26_plot_null_dat$condition, levels(rcp26_plot_null_dat$condition)[c(6,3,4,1,5,2)])

obs_vals <- rcp26_plot_null_dat %>% filter(condition==c('observed_6_15','observed_12_25','observed_24_45')) %>%
		group_by(condition) %>% summarise_all(mean)

rcp26_mod_null_perturbed <- ggplot(data=rcp26_plot_null_dat %>%
						filter(condition==c('random_phase_6_15','random_phase_12_25', 'random_phase_24_45')),
				aes(x=int_x, y=modularity, group=condition))+
		geom_violin(trim=T, show.legend=F, scale='width', aes(fill=condition), color='black')+
		geom_boxplot(width=0.2, show.legend=F, aes(fill=condition))+
		geom_point(data=obs_vals, aes(x=int_x, y=modularity, fill=condition), color='black', shape=23, size=5, show.legend=F)+
		scale_fill_manual(name = "", values = c("firebrick","dodgerblue3","darkorange2", 
						"firebrick","dodgerblue3","darkorange2"))+
		scale_color_manual(name = "", values = c("firebrick","dodgerblue3","darkorange2",  
						"firebrick","dodgerblue3","darkorange2"))+
		scale_x_continuous(name="", breaks=int.x,
				labels=c('RP\n6_15','OB\n6_15','RP\n12_25','OB\n12_25','RP\n24_45','OB\n24_45'),
				expand=c(0.05,0.05))+
		scale_y_continuous(name="Modularity", limits=c(0.88,0.98))+
		theme_bw()+
		theme(axis.text.x=element_text(angle=0, size=8),
				axis.text.y=element_text(size=10),
				axis.title.y=element_text(size=10),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank())

windows(width=4, height=5)
rcp26_mod_null_perturbed

#ggsave("~/MHWs/data/rcp26_perturb_modularity.pdf", useDingbats=FALSE, width = 4, height = 5)

pvals_mod_6_15 <- rcp26_plot_null_dat %>% filter(condition=='random_phase_6_15') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_6_15','modularity']>=
												as.numeric(obs_vals[obs_vals$condition=='observed_12_25','modularity'])))/1000))
pvals_mod_12_25 <- rcp26_plot_null_dat %>% filter(condition=='random_phase_12_25') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_12_25','modularity']>=
												as.numeric(obs_vals[obs_vals$condition=='observed_12_25','modularity'])))/1000))
pvals_mod_24_45 <- rcp26_plot_null_dat %>% filter(condition=='random_phase_24_45') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_24_45','modularity']>=
												as.numeric(obs_vals[obs_vals$condition=='observed_24_45','modularity'])))/1000))

rcp26_deg_null_perturbed <- ggplot(data=rcp26_plot_null_dat %>%
						filter(condition==c('random_phase_6_15','random_phase_12_25', 'random_phase_24_45')),
				aes(x=int_x, y=mean_deg, group=condition))+
		geom_violin(trim=T, show.legend=F, scale='width', aes(fill=condition), color='black')+
		geom_boxplot(width=0.2, show.legend=F, aes(fill=condition))+
		geom_point(data=obs_vals, aes(x=int_x, y=mean_deg, fill=condition), color='black', shape=23, size=5, show.legend=F)+
		scale_fill_manual(name = "", values = c("firebrick","dodgerblue3","darkorange2", 
						"firebrick","dodgerblue3","darkorange2"))+
		scale_color_manual(name = "", values = c("firebrick","dodgerblue3","darkorange2",  
						"firebrick","dodgerblue3","darkorange2"))+
		scale_x_continuous(name="", breaks=int.x,
				labels=c('RP\n6_15','OB\n6_15','RP\n12_25','OB\n12_25','RP\n24_45','OB\n24_45'),
				expand=c(0.05,0.05))+
		scale_y_continuous(name="Node degree")+
		theme_bw()+
		theme(axis.text.x=element_text(angle=0, size=8),
				axis.text.y=element_text(size=10),
				axis.title.y=element_text(size=10),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank())

windows(width=4, height=5)
rcp26_deg_null_perturbed

#ggsave("~/MHWs/data/rcp26_perturb_degree.pdf", useDingbats=FALSE, width = 4, height = 5)

pvals_deg_6_15 <- rcp26_plot_null_dat %>% filter(condition=='random_phase_6_15') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_6_15','mean_deg']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_6_15','mean_deg'])))/1000))
pvals_deg_12_25 <- rcp26_plot_null_dat %>% filter(condition=='random_phase_12_25') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_12_25','mean_deg']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_12_25','mean_deg'])))/1000))
pvals_deg_24_45 <- rcp26_plot_null_dat %>% filter(condition=='random_phase_24_45') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_24_45','mean_deg']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_24_45','mean_deg'])))/1000))

# ----- RCP 8.5 
load('~/MHWs/data/tda_rcp85_null_stats.RData')
load('~/MHWs/data/rcp85_optim1225_null_stats.RData')
load('~/MHWs/data/rcp85_optim3055_null_stats.RData')

rcp85_null_stats <- tda_rcp85_null_stats
rcp85_obs_dat <- rcp85_null_stats[[1]][1:1000,]
rcp85_phase_dat <- rcp85_null_stats[[2]][1:1000,]

rcp85_obs_optim1225 <- rcp85_optim1225_null_stats[[1]][1:1000,]
rcp85_phase_optim1225 <- rcp85_optim1225_null_stats[[2]][1:1000,]

rcp85_obs_optim3055 <- rcp85_optim3055_null_stats[[1]][1:1000,]
rcp85_phase_optim3055 <- rcp85_optim3055_null_stats[[2]][1:1000,]

rcp85_cond <- c(rep('random_phase_12_25', nrow(rcp85_phase_optim1225)),rep('observed_12_25', nrow(rcp85_phase_optim1225)),
		rep('random_phase_24_45', nrow(rcp85_phase_dat)), rep('observed_24_45', nrow(rcp85_obs_dat)),
		rep('random_phase_30_55', nrow(rcp85_phase_optim3055)), rep('observed_30_55', nrow(rcp85_phase_optim3055)))

int.x=c(0,1,3,4,6,7)
rcp85_plot_null_dat <- data.frame(int_x=c(rep(c(0,1), each=1000), rep(c(3,4), each=1000), rep(c(6,7), each=1000)),
		condition=rcp85_cond, rbind(rcp85_phase_optim1225,rcp85_obs_optim1225,
				rcp85_phase_dat,rcp85_obs_dat, rcp85_phase_optim3055, rcp85_obs_optim3055))

rcp85_plot_null_dat$condition <- factor(rcp85_plot_null_dat$condition)
rcp85_plot_null_dat$condition <- factor(rcp85_plot_null_dat$condition, levels(rcp85_plot_null_dat$condition)[c(6,3,4,1,5,2)])

obs_vals <- rcp85_plot_null_dat %>% filter(condition==c('observed_12_25','observed_24_45','observed_30_55')) %>%
		group_by(condition) %>% summarise_all(mean)

rcp85_mod_null_perturbed <- ggplot(data=rcp85_plot_null_dat %>%
						filter(condition==c('random_phase_12_25','random_phase_24_45', 'random_phase_30_55')),
				aes(x=int_x, y=modularity, group=condition))+
		geom_violin(trim=T, show.legend=F, scale='width', aes(fill=condition), color='black')+
		geom_boxplot(width=0.2, show.legend=F, aes(fill=condition))+
		geom_point(data=obs_vals, aes(x=int_x, y=modularity, fill=condition), color='black', shape=23, size=5, show.legend=F)+
		scale_fill_manual(name = "", values = c("darkorange2", "firebrick", "dodgerblue3",
						"darkorange2","firebrick", "dodgerblue3"))+
		scale_color_manual(name = "", values = c("darkorange2", "firebrick", "dodgerblue3",
						"darkorange2","firebrick", "dodgerblue3"))+
		scale_x_continuous(name="", breaks=int.x,
				labels=c('RP\n12_25','OB\n12_25','RP\n24_45','OB\n24_45','RP\n30_55','OB\n30_55'),
				expand=c(0.05,0.05))+
		scale_y_continuous(name="Modularity", limits=c(0.89,0.96))+
		theme_bw()+
		theme(axis.text.x=element_text(angle=0, size=8),
				axis.text.y=element_text(size=10),
				axis.title.y=element_text(size=10),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank())

windows(width=4, height=5)
rcp85_mod_null_perturbed

#ggsave("~/MHWs/data/rcp85_perturb_modularity.pdf", useDingbats=FALSE, width = 4, height = 5)

pvals_mod_12_25 <- rcp85_plot_null_dat %>% filter(condition=='random_phase_12_25') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_12_25','modularity']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_12_25','modularity'])))/1000))
pvals_mod_24_45 <- rcp85_plot_null_dat %>% filter(condition=='random_phase_24_45') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_24_45','modularity']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_24_45','modularity'])))/1000))
pvals_mod_30_55 <- rcp85_plot_null_dat %>% filter(condition=='random_phase_30_55') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_30_55','modularity']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_30_55','modularity'])))/1000))

rcp85_deg_null_perturbed <- ggplot(data=rcp85_plot_null_dat %>%
						filter(condition==c('random_phase_12_25','random_phase_24_45', 'random_phase_30_55')),
				aes(x=int_x, y=mean_deg, group=condition))+
		geom_violin(trim=T, show.legend=F, scale='width', aes(fill=condition), color='black')+
		geom_boxplot(width=0.2, show.legend=F, aes(fill=condition))+
		geom_point(data=obs_vals, aes(x=int_x, y=mean_deg, fill=condition), color='black', shape=23, size=5, show.legend=F)+
		scale_fill_manual(name = "", values = c("darkorange2", "firebrick", "dodgerblue3",
						"darkorange2","firebrick", "dodgerblue3"))+
		scale_color_manual(name = "", values = c("darkorange2", "firebrick", "dodgerblue3",
						"darkorange2","firebrick", "dodgerblue3"))+
		scale_x_continuous(name="", breaks=int.x,
				labels=c('RP\n12_25','OB\n12_25','RP\n24_45','OB\n24_45','RP\n30_55','OB\n30_55'),
				expand=c(0.05,0.05))+
		scale_y_continuous(name="Node degree", limits=c(2.3,3.3))+
		theme_bw()+
		theme(axis.text.x=element_text(angle=0, size=8),
				axis.text.y=element_text(size=10),
				axis.title.y=element_text(size=10),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank())

windows(width=4, height=5)
rcp85_deg_null_perturbed

#ggsave("~~/MHWs/data/rcp85_perturb_degree.pdf", useDingbats=FALSE, width = 4, height = 5)

pvals_deg_12_25 <- rcp85_plot_null_dat %>% filter(condition=='random_phase_12_25') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_12_25','mean_deg']<=
												as.numeric(obs_vals[obs_vals$condition=='observed_12_25','mean_deg'])))/1000))
pvals_deg_24_45 <- rcp85_plot_null_dat %>% filter(condition=='random_phase_24_45') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_24_45','mean_deg']>=
												as.numeric(obs_vals[obs_vals$condition=='observed_24_45','mean_deg'])))/1000))
pvals_deg_30_55 <- rcp85_plot_null_dat %>% filter(condition=='random_phase_30_55') %>%
		do(data.frame(p_value=length(which(.[.$condition=='random_phase_30_55','mean_deg']>=
												as.numeric(obs_vals[obs_vals$condition=='observed_30_55','mean_deg'])))/1000))

