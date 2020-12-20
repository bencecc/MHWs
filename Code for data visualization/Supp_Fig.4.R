#### ---- Plot the individual components of Supplementary Fig. 4, including ---- #### 
#### ---- the temporal connectivity matrices and node degree  ------------------ ####
#### ---- for Global (observed) MHWs and the RCP 2.6 scenario. ----------------- ####

require(reshape2)
require(ggplot2)
require(ggsci)
require(cmocean)
require(dplyr)
require(changepoint)
require(segmented)

#### ---- Global ---- #

# ---- Load the  data: 

# temporal connectivity matrix (data table) reduced by a factor 10 for plotting
load('~/MHWs/data/tcm_glob_red.RData') 
# list of results from random phase null models
load('~/MHWs/data/tda_glob_null_stats.RData') 
# days x # of years - useful for plotting
load('~/MHWs/data/glob_years_length.RData') 
#list with first object the full TCM data table and second object in the list is the normalized degree
# the full data table is heavy for plotting; use tcm_glob_red instead
load('~/MHWs/data/glob_res_mat.RData') 

#tcm_glob <- glob_res_mat[[1]] use this for plotting the original TCM, but it is a large matrix.
#tcm_glob_red provides a reduced matrix by factor 10 for visualization
norm_deg_glob <- glob_res_mat[[2]] #get node degree

plot_years <- c(1985,1990,1995,2000,2005,2010,2015)
int.time <- sapply(plot_years, function(x) {
			sel <- which(glob_years_length$year%in%x)
			sum(glob_years_length[1:sel,'no_days'])	
		})

melt_dat <- melt(as.matrix(tcm_glob_red))

tcm_glob <- ggplot(melt_dat, aes(x = Var1, y = rev(as.numeric(Var2)))) +
		geom_raster(aes(fill=value)) + coord_fixed() +
		scale_fill_gradientn(colours=cmocean('balance')(256))+
		scale_x_continuous(name="", breaks=c(int.time/10), expand=c(0,0), 
				labels=c('1985','1990','1995','2000','2005','2010','2015'))+
		scale_y_continuous(breaks=c(rev(int.time/10-min(int.time/10)+(1351.4-max(int.time/10)))), expand=c(0,0), 
				labels=c('1985','1990','1995','2000','2005','2010','2015'))+
		theme_bw(base_size = 12)+
		theme(axis.title = element_blank(),
				axis.text.x=element_text(size=10),
				axis.text.y=element_text(size=10),
				axis.line=element_blank(),
				panel.border=element_blank(),
				panel.grid=element_line(colour=NA))

windows(height=7, width=7)
plot(tcm_glob)

#ggsave("~/MHWs/data/tcm_glob.pdf", useDingbats=FALSE, width = 7, height = 7)

#### ---- plot node degree ---- ####
plot_years <- c(1982,1985,1990,1995,2000,2005,2010,2015,2018)
int.time <- sapply(plot_years, function(x) {
			sel <- which(glob_years_length$year%in%x)
			sum(glob_years_length[1:sel,'no_days'])	
		})

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

deg_glob_dat <- data.frame(x=1:length(norm_deg_glob), y=norm_deg_glob)

degree_glob_p <- ggplot(deg_null, aes(x=time, y=mean_deg))+
		scale_x_continuous(name="", breaks=c(int.time[2:8]), expand=c(0.05,0.05), 
				labels=c('1985','1990','1995','2000','2005','2010','2015'))+
		scale_y_continuous(name="Normalized degree", limits=c(-0.05,1.01), expand=c(0,0))+
		geom_line(data=deg_glob_dat, aes(x=x, y=y),
				col=pal_igv("alternating")(2)[1], alpha=0.5, linetype=4)+
		geom_ribbon(aes(ymin=ifelse(mean_deg-2*sd>=0,(mean_deg-2*sd),0),
						ymax=mean_deg+2*sd), fill=pal_igv("alternating")(2)[2], alpha=0.5)+
		geom_line(colour='white')+
		stat_smooth(data=deg_glob_dat, aes(x=x, y=y),
				col="#FC766AFF", se=F)+
		theme_bw(base_size = 12)+
		theme(panel.grid.major=element_line(colour=NA),
				axis.text.x=element_text(size=11),
				axis.text.y=element_text(size=11),
				axis.title.y=element_text(size=11),
				panel.grid.minor=element_line(colour=NA))


windows(width=5, height=3)
degree_glob_p

#ggsave("~/MHWs/data/node_glob.pdf", useDingbats=FALSE, width = 5, height = 3)

#### ---- RCP2.6 ---- #####
# temporal connectivity matrix (data table) reduced by a factor 10 for plotting
load('~/MHWs/data/tcm_rcp26_red.RData') 
# list of results from random phase null models
load('~/MHWs/data/tda_rcp26_null_stats.RData') 
# days x # of years - useful for plotting
load('~/MHWs/data/rcp26_years_length.RData') 
#list with first object the full TCM data table and second object in the list is the normalized degree
# the full data table is heavy for plotting; use tcm_rcp26_red instead
load('~/MHWs/data/rcp26_res_mat.RData') 

# normalized node degree
norm_deg_rcp26 <- rcp26_res_mat[[2]]

plot_years <- c(2010,2020,2030,2040,2050,2060,2070,2080,2090,2100)
int.time <- sapply(plot_years, function(x) {
			sel <- which(rcp26_years_length$year%in%x)
			sum(rcp26_years_length[1:sel,'no_days'])	
		})

#### ---- PLOTS ---- ####

melt_dat <- melt(as.matrix(tcm_rcp26_red))

tcm_rcp26_plot <- ggplot(melt_dat, aes(x = Var1, y = rev(as.numeric(Var2)), fill=value)) +
		geom_raster() + coord_fixed() +
		scale_x_continuous(breaks=int.time/10, expand=c(0,0), 
				labels=c('2010','2020','2030','2040','2050','2060','2070',
						'2080','2090','2100'))+
		scale_y_continuous(breaks=rev(int.time/10-min(int.time/10)+1), expand=c(0,0), 
				labels=(c('2010','2020','2030','2040','2050','2060','2070',
									'2080','2090','2100')))+
		scale_fill_gradientn(colours=cmocean('balance')(256))+
		theme_bw(base_size = 12)+
		theme(axis.title = element_blank(),
				axis.text.x=element_text(size=10),
				axis.text.y=element_text(size=10),
				axis.line=element_blank(),
				panel.border=element_blank(),
				panel.grid=element_line(colour=NA))

windows(height=7, width=7)
tcm_rcp26_plot
#ggsave("~/MHWs/data/tcm_rcp26.pdf", useDingbats=FALSE, width = 7, height = 7)

#### ---- Plot node degree RCP2.6 ---- ####
# ------- degree data from ensemble model
rcp26_df <- as.data.frame(tda_rcp26_null_stats[[4]])
rcp26_df_long <- stack(rcp26_df)
rcp26_deg_df <- data.frame(time=rep(1:34675, nrow(rcp26_df)), model=rep(1:nrow(rcp26_df), each=34675), degree=rcp26_df_long$value) 

deg_tmp <- apply(rcp26_df, 2, function(x) {
			mean=mean(x)
			sd=sd(x)
			c(mean, sd)
		})

deg_tmp <- t(deg_tmp)
deg_ind_rcp26 <- data.frame(time=1:nrow(deg_tmp), mean_deg=deg_tmp[,1],
		sd=deg_tmp[,2])

deg_dat_rcp26 <- data.frame(x=1:length(norm_deg_rcp26), y=norm_deg_rcp26)

rcp26_deg_plot <- ggplot(deg_ind_rcp26, aes(x=time, y=mean_deg))+
		scale_x_continuous(name="", breaks=c(int.time), expand=c(0.05,0.05), 
				labels=c('2010','2020','2030','2040','2050','2060','2070',
						'2080','2090','2100'))+
		scale_y_continuous(name="Normalized degree", limits=c(-0.05,1.01), expand=c(0,0))+
		geom_line(data=deg_dat_rcp26, aes(x=x, y=y),
				col=pal_igv("alternating")(2)[1], alpha=0.3, linetype=4)+
		geom_ribbon(aes(ymin=ifelse(mean_deg-2*sd>=0,(mean_deg-2*sd),0),
						ymax=mean_deg+2*sd), fill=pal_igv("alternating")(2)[2], alpha=0.5)+
		geom_line(colour='white')+
		stat_smooth(data=deg_dat_rcp26, aes(x=x, y=norm_deg_rcp26),
				col="#FC766AFF", se=F)+
		theme_bw(base_size = 12)+
		theme(panel.grid.major=element_line(colour=NA),
				axis.text.x=element_text(size=11),
				axis.text.y=element_text(size=11),
				axis.title.y=element_text(size=11),
				panel.grid.minor=element_line(colour=NA))


windows(height=3, width=5)
rcp26_deg_plot

#ggsave("~/MHWs/data/node_rcp26.pdf", useDingbats=FALSE, width = 5, height = 3)



