#### ---- Plot the individual components of Fig. 3, including ---- #### 
#### ---- the temporal connectivity matrices and node degree  ---- ####
#### ---- for Historical and RCP 8.5 scenarios. ------------------ ####

###############################################################################
require(reshape2)
require(ggplot2)
require(mgcv)
require(ggsci)
require(cmocean)
require(dplyr)
require(changepoint)
require(segmented)
require(boot)

#### ---- HISTORICAL ---- ####

# ---- Load the data: 
# temporal connectivity matrix (data table) reduced by a factor 10 for plotting
load('~/MHWs/data/tcm_hist_red.RData') 
#list with first object the full TCM data table and second object in the list is the normalized degree
# the full data table is heavy for plotting; use tcm_hist_red instead
load('~/MHWs/data/hist_res_mat.RData') 
# list of results from random phase null models
load('~/MHWs/data/tda_hist_null_stats.RData') 
# load node degree from individual hist models with mean and sd from null distirbutions
# and collapse year 
load('~/MHWs/data/hist_degree_ind.RData')
# load mean and sd null model data for CIs of node degree
load('~/MHWs/data/hist_deg_null.RData')


# normalized node degree
norm_deg_hist <- hist_res_mat[[2]] 

# data in long format to plot the TCM
melt_dat <- melt(as.matrix(tcm_hist_red))

tmp=seq(as.Date("1861/1/1"), as.Date("2005/12/31"), by = "day")
years_all <- data.frame(year=format(tmp, "%Y"))
hist_years_length <- years_all %>% group_by(year) %>% summarise(no_days=n())

plot_years <- c(1880,1900,1920,1940,1960,1980,2000)
int.time <- sapply(plot_years, function(x) {
			sel <- which(hist_years_length$year%in%x)
			sum(hist_years_length[1:sel,'no_days'])	
		})


# ---- detect change points -> time scales of change in the TCM 
change_point_hist <- cpt.mean(norm_deg_hist, method='PELT', penalty = "Manual", pen.value = "7.5 * log(n)")
segms <- cpts(change_point_hist) #segments of interest for plotting 1, 3, 5

# ----------------------------------------------------------------------------- #

# ---- determine the temporal scales of coherence: 63, 26 and 27 yrs for the three major blocks, respectively
vec <- c(1,segms,52960)
diff(vec)/365
# determine timescales of coherence for the voids in the TCM: from 1 day to two years
load('spread_yr_hist.RData')
range(c(unique(spread_yr_hist[c(segms[1]:segms[2]-12),1])[1:70],
unique(spread_yr_hist[c(segms[3]:segms[4]),1])[1:64],
unique(spread_yr_hist[c(segms[5]:nrow(spread_yr_hist)),1])))

# --------------------------------------------------- #

# ---- determine the years of coherence
# first block (sync)
range(as.vector(years_all[1:segms[1],]))
# second block (void)
range(as.vector(years_all[segms[1]:segms[2],]))
# third block (sync)
range(as.vector(years_all[segms[2]:segms[3],]))
# fourth block (void)
range(as.vector(years_all[segms[3]:segms[4],]))
# fifth block (sync)
range(as.vector(years_all[segms[4]:segms[5],]))
# sixth block (void)
range(as.vector(years_all[segms[5]:nrow(years_all),]))

# --------------------------------------------- #

# build polygons to represent the temporal scale of coherence identified by the changepoint procedure
s1 <- c(1,segms[1]/10) #divide by 10 for plotting
s2 <- c(segms[2]/10,segms[3]/10)
s3 <- c(segms[4]/10,segms[5]/10) #timescales correspond to length of s1-s3
top11 <- rev(as.numeric(melt_dat$Var2))[1]
top12 <- top11-s1[2]
top21 <- top12-(s2[1]-s1[2])
top22 <- top21-(s2[2]-s2[1])
top31 <- top22-(s3[1]-s2[2])
top32 <- top31-(s3[2]-s3[1])

p1 <- data.frame(x=c(s1[1]+20, s1[2], s1[2], s1[1]+20, s1[1]+20),
		y=c(top11-20, top11-20, top12, top12, top11-20))
p2 <- data.frame(x=c(s2[1], s2[2], s2[2], s2[1], s2[1]),
		y=c(top21, top21, top22, top22, top21))
p3 <- data.frame(x=c(s3[1], s3[2], s3[2], s3[1], s3[1]),
		y=c(top31, top31, top32, top32, top31))

# data in long format to plot the TCM
melt_dat <- melt(as.matrix(tcm_hist_red))

tcm_hist <- ggplot(melt_dat, aes(x = Var1, y = rev(as.numeric(Var2)))) +
		geom_raster(aes(fill=value)) + coord_fixed() +
		scale_fill_gradientn(colours=cmocean('balance')(256))+
		geom_polygon(data=p1, aes(x=x,y=y), col='yellow', fill=NA, size=0.5)+
		geom_polygon(data=p2, aes(x=x,y=y), col='yellow', fill=NA, size=0.5)+
		geom_polygon(data=p3, aes(x=x,y=y), col='yellow', fill=NA, size=0.5)+
		scale_x_continuous(name="", breaks=c(1,int.time/10), expand=c(0,0), 
				labels=c('1861','1880','1900','1920','1940','1960','1980','2000'))+
		scale_y_continuous(breaks=c(5296,rev(int.time/10-min(int.time/10)+(5296-max(int.time/10)))), expand=c(0,0), 
				labels=c('1861','1880','1900','1920','1940','1960','1980','2000'))+
		theme_bw(base_size = 12)+
		theme(axis.title = element_blank(),
				axis.text.x=element_text(size=10),
				axis.text.y=element_text(size=10),
				axis.line=element_blank(),
				panel.border=element_blank(),
				panel.grid=element_line(colour=NA))

windows(height=7, width=7)
plot(tcm_hist)

#ggsave("~/MHWs/data/tcm_hist.pdf", useDingbats=FALSE, width = 7, height = 7)

##### --- PLOT NODE DEGREE ---- #####
# ------- degree data from ensemble model
norm_deg_hist <- hist_res_mat[[2]] # normalized node degree
hist_deg_dat <- data.frame(x=as.numeric(names(norm_deg_hist)), y=norm_deg_hist)

# ------------ Identify year when observed degree is undistinguisheable from null
deg_vec <- predict(gam(y~s(x, bs = "cs"), data=hist_deg_dat))
deg_diff <- as.numeric(deg_vec)-(hist_deg_null$mean_deg+2*hist_deg_null$sd)
min_year_deg <- which(deg_diff>0)
shift_time <- min_year_deg[length(min_year_deg)] + 1
# identify year of collapse
collapse_year_hist <- years_all %>% slice(shift_time)

# ------- normalized node degree data from individual null models (use 500 randomizations)
hist_df <- as.data.frame(tda_hist_null_stats[[4]][1:500,])
hist_df_long <- stack(hist_df)
hist_deg_df <- data.frame(time=rep(1:52960, nrow(hist_df)), model=rep(1:nrow(hist_df), each=52960), degree=hist_df_long$value) 

deg_tmp <- apply(hist_df, 2, function(x) {
			mean=mean(x)
			sd=sd(x)
			c(mean, sd)
		})

deg_tmp <- t(deg_tmp)
hist_deg_null <- data.frame(time=1:nrow(deg_tmp), mean_deg=deg_tmp[,1],
		sd=deg_tmp[,2])

# --------- Get the median and bootstrapped standard error of the year of collapse
# --------- from individual models
coll_yr_hist <- hist_degree_ind %>% group_by(model) %>% select(collapse_year) %>% distinct()
shift_median_hist <- median(coll_yr_hist$collapse_year, na.rm=T)
median_shift_year_hist <- years_all %>% slice(shift_median_hist) 
med_boot_hist <- boot(data = na.exclude(coll_yr_hist$collapse_year),
		statistic = function(x,i) median(x[i]),R = 10000)
boot_se_hist <- sd(med_boot_hist$t[,1])

degree_hist_p <- ggplot(hist_deg_null, aes(x=time, y=mean_deg))+
		scale_x_continuous(name="", breaks=c(1,int.time), expand=c(0.05,0.05), 
				labels=c('1861','1880','1900','1920','1940','1960','1980','2000'))+
		scale_y_continuous(name="Normalized degree", limits=c(-0.05,1.01), expand=c(0,0))+
		geom_line(data=hist_deg_dat, aes(x=x, y=y),
				col=pal_igv("alternating")(2)[1], alpha=0.3, linetype=4)+
		geom_rect(aes(xmin=shift_median_hist-boot_se_hist,
						xmax=ifelse((shift_median_hist+boot_se_hist)>nrow(hist_deg_null),
								nrow(hist_deg_null),(shift_median_hist+boot_se_hist)),
						ymin=-0.05, ymax=0), fill="#FC766AFF", alpha=1)+
		geom_ribbon(aes(ymin=ifelse(mean_deg-2*sd>=0,(mean_deg-2*sd),0),
						ymax=mean_deg+2*sd), fill=pal_igv("alternating")(2)[2], alpha=0.5)+
		geom_line(colour='white')+
		scale_color_gradientn(colours=cmocean('matter')(256))+
		stat_smooth(data=hist_deg_dat, aes(x=x, y=norm_deg_hist), method="gam", formula = y ~ s(x, bs = "cs"),
				col="#FC766AFF", se=F)+
		geom_segment(x=shift_time, xend=shift_time, y=1.05, yend=-0.05, col="grey40", size=0.7)+
		geom_segment(x=shift_median_hist, xend=shift_median_hist, y=0, yend=-0.05, col="green", size=1)+
		labs(x="", y="Normalized degree")+
		theme_bw(base_size = 12)+
		theme(panel.grid.major=element_line(colour=NA),
				axis.text.x=element_text(size=11),
				axis.text.y=element_text(size=11),
				axis.title.y=element_text(size=11),
				panel.grid.minor=element_line(colour=NA))

windows(width=5, height=3)
degree_hist_p

#ggsave("~/MHWs/data/node_hist.pdf", width = 5, height = 3)

#### ---- RCP85 ---- #####
# ---- Load the data: 
# temporal connectivity matrix (data table) reduced by a factor 10 for plotting
load('~/MHWs/data/tcm_rcp85_red.RData') 
#list with first object the full TCM data table and second object in the list is the normalized degree
# the full data table is heavy for plotting; use tcm_rcp85_red instead
load('~/MHWs/data/rcp85_res_mat.RData') 
# list of results from random phase null models
load('~/MHWs/data/tda_rcp85_null_stats.RData') 
# load node degree from individual rcp85 models with mean and sd from null distirbutions
# and collapse year 
load('~/MHWs/data/rcp85_degree_ind.RData')
# load mean and sd null model data for CIs of node degree
load('~/MHWs/data/rcp85_deg_null.RData')


# normalized node degree
norm_deg_rcp85 <- rcp85_res_mat[[2]]

tmp=seq(as.Date("2006/1/1"), as.Date("2100/12/31"), by = "day")
years_all <- data.frame(year=format(tmp, "%Y"))
rcp85_years_length <- years_all %>% group_by(year) %>% summarise(no_days=n())

plot_years <- c(2010,2020,2030,2040,2050,2060,2070,2080,2090,2100)
int.time <- sapply(plot_years, function(x) {
			sel <- which(rcp85_years_length$year%in%x)
			sum(rcp85_years_length[1:sel,'no_days'])	
		})

melt_dat <- melt(as.matrix(tcm_rcp85_red))
change_point_rcp85 <- cpt.mean(norm_deg_rcp85, method="PELT", penalty="Manual", pen.value="2.5*log(n)") 
segms_rcp85 <- cpts(change_point_rcp85)

s1_rcp85 <- c(segms_rcp85[1]/10, nrow(tcm_rcp85_red )) #divide by 10 for plotting
top11_rcp85 <- rev(as.numeric(melt_dat$Var2))[1]-s1_rcp85[1]
top12_rcp85 <- top11_rcp85 - (s1_rcp85[2]-s1_rcp85[1]) #

p1_rcp85 <- data.frame(x=c(s1_rcp85[1], s1_rcp85[2]-10, s1_rcp85[2]-10, s1_rcp85[1], s1_rcp85[1]),
		y=c(top11_rcp85, top11_rcp85, top12_rcp85+10, top12_rcp85+10, top11_rcp85))

tcm_rcp85_plot <- ggplot(melt_dat, aes(x = Var1, y = rev(as.numeric(Var2)), fill=value)) +
		geom_raster() + coord_fixed() +
#		geom_polygon(data=p1_rcp85, aes(x=x, y=y), fill=NA, col='yellow', size=0.5)+
		scale_fill_gradientn(colours=cmocean('balance')(256))+
		scale_x_continuous(breaks=int.time/10, expand=c(0,0), 
				labels=c('2010','2020','2030','2040','2050','2060','2070',
						'2080','2090','2100'))+
		scale_y_continuous(breaks=rev(int.time/10-min(int.time/10)+1), expand=c(0,0), 
				labels=(c('2010','2020','2030','2040','2050','2060','2070',
						'2080','2090','2100')))+
		theme_bw(base_size = 12)+
		theme(axis.title = element_blank(),
				axis.text.x=element_text(size=10),
				axis.text.y=element_text(size=10),
				axis.line=element_blank(),
				panel.border=element_blank(),
				panel.grid=element_line(colour=NA))

windows(height=7, width=7)
tcm_rcp85_plot

#ggsave("~/MHWs/data/tcm_rcp85.pdf", useDingbats=FALSE, width = 7, height = 7)


#### ---- Plot Node Degree ---- ####
# ------- degree data from ensemble model
norm_deg_rcp85 <- rcp85_res_mat[[2]]
rcp85_deg_dat <- data.frame(x=as.numeric(names(norm_deg_rcp85)), y=norm_deg_rcp85)

# ------------------- identify year when observed degree is undistinguisheable from null ####
deg_vec <- predict(gam(y~s(x, bs = "cs"), data=rcp85_deg_dat))
deg_diff <- as.numeric(deg_vec)-(rcp85_deg_null$mean_deg+2*rcp85_deg_null$sd)
#min_year_deg <- which(deg_diff==min(deg_diff[deg_diff>0]))
min_year_deg <- which(deg_diff>0)
shift_time_rcp85 <- min_year_deg[1]
collapse_year_rcp85 <- years_all %>% slice(shift_time_rcp85)

# ------- normalized node degree data from individual null models (use 500 randomizations)
rcp85_df <- as.data.frame(tda_rcp85_null_stats[[4]][1:500,])
rcp85_df_long <- stack(rcp85_df)
rcp85_deg_df <- data.frame(time=rep(1:34675, nrow(rcp85_df)), model=rep(1:nrow(rcp85_df), each=34675), degree=rcp85_df_long$value) 

deg_tmp <- apply(rcp85_df, 2, function(x) {
			mean=mean(x)
			sd=sd(x)
			c(mean, sd)
		})

deg_tmp <- t(deg_tmp)
rcp85_deg_null <- data.frame(time=1:nrow(deg_tmp), mean_deg=deg_tmp[,1],
		sd=deg_tmp[,2])

# --------- Get the median and bootstrapped standard error of the year of collapse
# --------- from individual models
coll_yr_rcp85 <- rcp85_degree_ind %>% group_by(model) %>% select(collapse_year) %>% distinct()
shift_median_rcp85 <- median(coll_yr_rcp85$collapse_year, na.rm=T)
median_shift_year_rcp85 <- years_all %>% slice(shift_median_rcp85) 

med_boot_rcp85 <- boot(data = na.exclude(coll_yr_rcp85$collapse_year),
		statistic = function(x,i) median(x[i]),R = 10000)
boot_se_rcp85 <- sd(med_boot_rcp85$t[,1])

# ------------ Plot 
degree_rcp85_p <- ggplot(rcp85_deg_null, aes(x=time, y=mean_deg))+
		scale_x_continuous(name="", breaks=int.time, expand=c(0.05,0.05), 
				labels=c('2010','2020','2030','2040','2050','2060','2070',
						'2080','2090','2100'))+	
		scale_y_continuous(name="Normalized degree", limits=c(-0.05,1.01), expand=c(0,0))+
		geom_line(data=rcp85_deg_dat, aes(x=x, y=y),
				col=pal_igv("alternating")(2)[1], alpha=0.3, linetype=4)+
		geom_rect(aes(xmin=shift_median_rcp85-boot_se_rcp85,
						xmax=ifelse((shift_median_rcp85+boot_se_rcp85)>nrow(rcp85_deg_null),
								nrow(rcp85_deg_null),(shift_median_rcp85+boot_se_rcp85)),
						ymin=-0.05, ymax=0), fill="#FC766AFF", alpha=1)+
		geom_ribbon(aes(ymin=ifelse(mean_deg-2*sd>=0,(mean_deg-2*sd),0),
						ymax=mean_deg+2*sd), fill=pal_igv("alternating")(2)[2], alpha=0.5)+
		geom_line(colour='white')+
		scale_color_gradientn(colours=cmocean('matter')(256))+
		stat_smooth(data=rcp85_deg_dat, aes(x=x, y=norm_deg_rcp85), method="gam", formula = y ~ s(x, bs = "cs"),
				col="#FC766AFF", se=F)+
		geom_segment(x=shift_time_rcp85, xend=shift_time_rcp85, y=1.05, yend=-0.05, col="grey40", size=0.7)+
		geom_segment(x=shift_median_rcp85, xend=shift_median_rcp85, y=0, yend=-0.05, col="green", size=1)+
		labs(x="", y="Normalized degree")+
		theme_bw(base_size = 12)+
		theme(panel.grid.major=element_line(colour=NA),
				axis.text.x=element_text(size=11),
				axis.text.y=element_text(size=11),
				axis.title.y=element_text(size=11),
				panel.grid.minor=element_line(colour=NA))

windows(width=5, height=3)
degree_rcp85_p

#ggsave("~/MHWs/data/node_rcp85.pdf", width = 5, height = 3)

