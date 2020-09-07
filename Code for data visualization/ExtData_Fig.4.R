# ---------- Plot node degree of temporal connectivity ----------------- #
# ---------- of individual CMIP5 models (Historical and RCP 8.5) ------- #
# ---------- matrices with random phase null models -------------------- #

require(ggplot2)
require(ggsci)
require(cmocean)
require(dplyr)
require(boot)

# -------- Historical ---------------#
# load node degree from individual hist models with mean and sd from null distributions
# and collapse year 
load('~/MHWs/data/hist_degree_ind.RData')

tmp=seq(as.Date("1861/1/1"), as.Date("2005/12/31"), by = "day")
years_all <- data.frame(year=format(tmp, "%Y"))
hist_years_length <- years_all %>% group_by(year) %>% summarise(no_days=n())

plot_years <- c(1880,1900,1920,1940,1960,1980,2000)
int.time <- sapply(plot_years, function(x) {
			sel <- which(hist_years_length$year%in%x)
			sum(hist_years_length[1:sel,'no_days'])	
		})

coll_yr_hist <- hist_degree_ind %>% group_by(model) %>% select(collapse_year) %>% distinct()
shift_median_hist <- median(coll_yr_hist$collapse_year, na.rm=T)
median_shift_year_hist <- years_all %>% slice(shift_median_hist) 

med_boot_hist <- boot(data = na.exclude(coll_yr_hist$collapse_year),
		statistic = function(x,i) median(x[i]),R = 10000)
boot_se_hist <- sd(med_boot_hist$t[,1])

hist_ind_p1 <- ggplot(hist_degree_ind, aes(x=time, y=degree))+
		scale_x_continuous(name="", breaks=int.time, expand=c(0.05,0.05), 
				labels=c('1880','1900','1920','1940','1960','1980','2000'))+
		scale_y_continuous(name="Normalized degree", limits=c(-0.05,1.01), expand=c(0,0))+
		geom_line(col=pal_igv("alternating")(2)[1], alpha=0.3, linetype=4)+
		geom_ribbon(aes(ymin=mean_null-2*sd_null, ymax=mean_null+2*sd_null), fill=pal_igv("alternating")(2)[2], alpha=0.5)+
		geom_line(aes(x=time, y=mean_null),colour='white')+
		stat_smooth(method="gam", formula = y ~ s(x, bs = "cs"),
				col="#FC766AFF", se=F)+
		scale_color_gradientn(colours=cmocean('matter')(256))+
		geom_segment(x=hist_degree_ind$collapse_year, xend=hist_degree_ind$collapse_year,
				y=1.05, yend=-0.05, col="grey40", size=0.7)+
		facet_wrap(~model, scales='fixed',ncol=3)+
		labs(x="", y="Normalized degree")+
		theme_bw(base_size = 12)+
		theme(panel.grid.major=element_line(colour=NA),
				axis.text.x=element_text(size=8),
				axis.text.y=element_text(size=8),
				axis.title.y=element_text(size=8),
				strip.text = element_text(size=8),
				strip.background = element_rect(size=0.5),
				panel.grid.minor=element_line(colour=NA))

windows(width=9, height=5)
print(hist_ind_p1)

#ggsave("~/MHWs/data/node_hist_ind.pdf", width = 9, height = 5)

# -------- RCP 8.5 data ------------ #

# load node degree from individual rcp85 models with mean and sd from null distributions
# and collapse year 
load('~/MHWs/data/rcp85_degree_ind.RData')

tmp=seq(as.Date("2006/1/1"), as.Date("2100/12/31"), by = "day")
years_all <- data.frame(year=format(tmp, "%Y"))
rcp85_years_length <- years_all %>% group_by(year) %>% summarise(no_days=n())

plot_years <- c(2010,2020,2030,2040,2050,2060,2070,2080,2090,2100)
int.time <- sapply(plot_years, function(x) {
			sel <- which(rcp85_years_length$year%in%x)
			sum(rcp85_years_length[1:sel,'no_days'])	
		})

coll_yr_rcp85 <- rcp85_degree_ind %>% group_by(model) %>% select(collapse_year) %>% distinct()
shift_median_rcp85 <- median(coll_yr_rcp85$collapse_year, na.rm=T)
median_shift_year_rcp85 <- years_all %>% slice(shift_median_rcp85) 
med_boot_rcp85 <- boot(data = na.exclude(coll_yr_rcp85$collapse_year),
		statistic = function(x,i) median(x[i]),R = 10000)
boot_se_rcp85 <- sd(med_boot_rcp85$t[,1])


ind_p1 <- ggplot(rcp85_degree_ind, aes(x=time, y=degree))+
		scale_x_continuous(name="", breaks=int.time, expand=c(0.05,0.05), 
				labels=c('2010','2020','2030','2040','2050','2060','2070',
						'2080','2090','2100'))+	
		scale_y_continuous(name="Normalized degree", limits=c(-0.05,1.01), expand=c(0,0))+
		geom_line(col=pal_igv("alternating")(2)[1], alpha=0.3, linetype=4)+
		geom_ribbon(aes(ymin=ifelse(mean_null-2*sd_null>=0,(mean_null-2*sd_null),0),
						ymax=mean_null+2*sd_null), fill=pal_igv("alternating")(2)[2], alpha=0.5)+
		geom_line(aes(x=time, y=mean_null),colour='white')+
		stat_smooth(method="gam", formula = y ~ s(x, bs = "cs"),
				col="#FC766AFF", se=F)+
		scale_color_gradientn(colours=cmocean('matter')(256))+
		geom_segment(x=rcp85_degree_ind$collapse_year, xend=rcp85_degree_ind$collapse_year,
				y=1.05, yend=-0.05, col="grey40", size=0.7)+
		facet_wrap(~model, scales='fixed',ncol=3)+	
		labs(x="", y="Normalized degree")+
		theme_bw(base_size = 12)+
		theme(panel.grid.major=element_line(colour=NA),
				axis.text.x=element_text(size=8),
				axis.text.y=element_text(size=8),
				axis.title.y=element_text(size=8),
				strip.text = element_text(size=8),
				strip.background = element_rect(size=0.5),
				panel.grid.minor=element_line(colour=NA))

windows(width=9, height=5)
ind_p1

#ggsave("~/MHWs/data/node_rcp85_ind.pdf", width = 9, height = 5)








