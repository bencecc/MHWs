#### ----------- Supplementary Fig. 9. Distance distribution ----- ####
#### ----------- of significantly synchronized connections. ------ ####

require(ggplot2)
require(scales)
require(dplyr)
require(tidyr)

load("~/MHWs/data/es_plotdat_tau10.RData")
load("~/MHWs/data/es_plotdat_tau30.RData")
load("~/MHWs/data/es_plotdat_tau10bis.RData")
load("~/MHWs/data/es_plotdat_tau30bis.RData")

#### ---- Plot components of Extended Data Fig. 10 ---- ####
#### ---- Supplementary figure: plot the analogue of Fig 4 of main manuscript with tau = 30 ---- ####
es_dist_dat_tau30 <- es_plotdat_tau30

dist.pos_tau30 <- c(seq(1, 200, by=10), seq(201, 315, by=20), seq(316, 380, by=35),
		seq(381,460, by=40), seq(461, 800, by=50), seq(801,1100, by=70), seq(1101,2000, by=120))
dist_dat_tau30 <- es_dist_dat_tau30[dist.pos_tau30,] %>% gather(reference, pdf, c(2:7))
dist_dat_tau30$line_type <- rep("s", nrow(dist_dat_tau30))

power_dat_tau30 <- es_dist_dat_tau30 %>% filter(distance>=500&distance<=3600) %>%
		select(distance, pdf=gcdist_declust30_sat_1982_2018,
				pdf1=gcdist_declust30_1891_1900,
				pdf2=gcdist_declust30_1991_2000,
				gcdist_all)
power_fit_tau30 <- lm(log10(pdf) ~ log10(distance), data=power_dat_tau30)
summary(power_fit_tau30)
x_range_tau30 <- seq(500, 10^4.29, length.out=2000)
pred_lines_tau30 <- data.frame(distance=x_range_tau30, pdf=10^(predict(power_fit_tau30, data.frame(distance=x_range_tau30))))

power_fit1_tau30 <- lm(log10(pdf1) ~ log10(distance), data=power_dat_tau30)
summary(power_fit1_tau30)

power_fit2_tau30 <- lm(log10(pdf2) ~ log10(distance), data=power_dat_tau30)
summary(power_fit2_tau30)

gc_diff_tau30 <- es_dist_dat_tau30 %>% filter(distance>10^3.4, distance<10^3.6) %>%
		mutate(diff=log10(gcdist_declust30_sat_1982_2018)-log10(gcdist_all)) %>%
		select(distance, diff, gcdist_all)
min.sel_tau30 <- which(gc_diff_tau30$diff<0)[1]
inters_dist_tau30 <- gc_diff_tau30[(min.sel_tau30)-1,"distance"]
inters_y_tau30 <- es_dist_dat_tau30[which(es_dist_dat_tau30[,"distance"]>inters_dist_tau30)[1]-1, 2]
y_line_tau30 <- data.frame(x=rep(inters_dist_tau30, 20), y=seq(0, inters_y_tau30, length.out=20))

gc_diff1_tau30 <- es_dist_dat_tau30 %>% filter(distance>10^3.6, distance<10^3.8) %>%
		mutate(diff=log10(gcdist_declust30_2081_2090)-log10(gcdist_all)) %>%
		select(distance, diff, gcdist_all)
min.sel1_tau30 <- which(gc_diff1_tau30$diff<0)[1]
inters_dist1_tau30 <- gc_diff1_tau30[(min.sel1_tau30)-1,"distance"]
inters_y1_tau30 <- es_dist_dat_tau30[which(es_dist_dat_tau30[,"distance"]>inters_dist1_tau30)[1]-1, 2]
y_line1_tau30 <- data.frame(x=rep(inters_dist1_tau30, 20), y=seq(0, inters_y1_tau30, length.out=20))

line_dat_tau30 <- dist_dat_tau30 %>% filter(reference!="gcdist_declust30_sat_1982_2018")
point_dat_tau30 <- dist_dat_tau30 %>% filter(reference=="gcdist_declust30_sat_1982_2018")
point_dat_tau30$aest <- "s"

p_dist_tau30 <- ggplot(data=line_dat_tau30, aes(x=distance, y=pdf))
pp_dist_tau30 <- p_dist_tau30 +
		geom_line(aes(x=distance, y=pdf, col=reference, size=reference))+
		geom_point(data=point_dat_tau30, aes(x=distance, y=pdf, col=aest), shape=1, size=3)+
		geom_line(data=y_line_tau30, aes(x=x, y=y), linetype="solid", color = "#CC79A7", size=0.5)+
		geom_line(data=y_line1_tau30, aes(x=x, y=y), linetype="solid", color = "#E69F00", size=0.5)+
		geom_line(data=pred_lines_tau30, aes(x=distance, y=pdf),
				col="grey40", size=0.5, linetype="dashed")+
		scale_size_manual(name=NULL, values=c(1, rep(0.8,4), 3), guide="none")+
		scale_color_manual(name=NULL, labels=c("Great-circle KDE", "Historical (1891-1900)", "Historical (1991-2000)",
						"RCP8.5 (2041-2050)", "RCP8.5 (2081-2090)","Satellite (1982-2018)"),
				values = c("black","#00AFBB","#0072B2","#999999","#D55E00","#CC6666"),
				guide = guide_legend(override.aes = list(linetype = c(rep("solid", 5), "blank"), shape = c(rep(NA,5), 1))))+
		scale_x_log10(name="Distance (km)", limits=c(1e2, 1e5), breaks = trans_breaks("log10", function(x) 10^x),
				labels = trans_format("log10", math_format(10^.x)))+
		scale_y_log10(limits=c(10^-6.2, 10^-3), name = "Probability density", breaks = trans_breaks("log10", function(x) 10^x),
				labels = trans_format("log10", math_format(10^.x)))+
		theme_classic()+
		theme(legend.key.size = unit(0.15, "cm"),
				legend.position=c(0.8,0.9),
				legend.spacing = unit(0.7, "cm"))

windows(height=4, width=6)
pp_dist_tau30

#ggsave("~/MHWs/data/event_sync_tau30.pdf", useDingbats=FALSE, width = 6, height = 4)

#### ---- Supplementary figure with alternative periods for ES (tau10) ---- ####
es_dist_dat_tau10bis <- es_plotdat_tau10bis

dist.pos_tau10bis <- c(seq(1, 200, by=20), seq(201, 315, by=40), seq(316, 380, by=80),
		seq(381,560, by=100), seq(561,1100, by=80), seq(1101,2000, by=120))
dist_dat_tau10bis <- es_dist_dat_tau10bis[dist.pos_tau10bis,] %>% gather(reference, pdf, c(2:7))
dist_dat_tau10bis$line_type <- rep("s", nrow(dist_dat_tau10bis))

# to fit a power law: select a range between 10^2.6021 and 10^3.5 corresponding to 400-3162 km
power_dat_tau10bis <- es_dist_dat_tau10bis %>% filter(distance>=500&distance<=3600) %>%
		select(distance, pdf=gcdist_declust_sat_1982_2018,
				pdf1=gcdist_declust_1861_1923,
				pdf2=gcdist_declust_1991_2005,
				gcdist_all)
power_fit_tau10bis <- lm(log10(pdf) ~ log10(distance), data=power_dat_tau10bis)
summary(power_fit_tau10bis)
x_range_tau10bis <- seq(500, 10^4.29, length.out=2000)
pred_lines_tau10bis <- data.frame(distance=x_range_tau10bis, pdf=10^(predict(power_fit_tau10bis, data.frame(distance=x_range_tau10bis))))

power_fit1_tau10bis <- lm(log10(pdf1) ~ log10(distance), data=power_dat_tau10bis)
summary(power_fit1_tau10bis)

power_fit2_tau10bis <- lm(log10(pdf2) ~ log10(distance), data=power_dat_tau10bis)
summary(power_fit2_tau10bis)

gc_diff_tau10bis <- es_dist_dat_tau10bis %>% filter(distance>10^3.4, distance<10^3.6) %>%
		mutate(diff=log10(gcdist_declust_sat_1982_2018)-log10(gcdist_all)) %>%
		select(distance, diff, gcdist_all)
min.sel_tau10bis <- which(gc_diff_tau10bis$diff<0)[1]
inters_dist_tau10bis <- gc_diff_tau10bis[(min.sel_tau10bis)-1,"distance"]
inters_y_tau10bis <- es_dist_dat_tau10bis[which(es_dist_dat_tau10bis[,"distance"]>inters_dist_tau10bis)[1]-1, 2]
y_line_tau10bis <- data.frame(x=rep(inters_dist_tau10bis, 20), y=seq(0, inters_y_tau10bis, length.out=20))

gc_diff1_tau10bis <- es_dist_dat_tau10bis %>% filter(distance>10^3.6, distance<10^4) %>%
		mutate(diff=log10(gcdist_declust_2071_2100)-log10(gcdist_all)) %>%
		select(distance, diff, gcdist_all)
min.sel1_tau10bis <- which(gc_diff1_tau10bis$diff<0)[1]
inters_dist1_tau10bis <- gc_diff1_tau10bis[(min.sel1_tau10bis)-1,"distance"]
inters_y1_tau10bis <- es_dist_dat_tau10bis[which(es_dist_dat_tau10bis[,"distance"]>inters_dist1_tau10bis)[1]-1, 2]
y_line1_tau10bis <- data.frame(x=rep(inters_dist1_tau10bis, 20), y=seq(0, inters_y1_tau10bis, length.out=20))

line_dat_tau10bis <- dist_dat_tau10bis %>% filter(reference!="gcdist_declust_sat_1982_2018")
point_dat_tau10bis <- dist_dat_tau10bis %>% filter(reference=="gcdist_declust_sat_1982_2018")
point_dat_tau10bis$aest <- "s"

p_dist_tau10bis <- ggplot(data=line_dat_tau10bis, aes(x=distance, y=pdf))
pp_dist_tau10bis <- p_dist_tau10bis +
		geom_line(aes(x=distance, y=pdf, col=reference, size=reference))+
		geom_point(data=point_dat_tau10bis, aes(x=distance, y=pdf, col=aest), shape=1, size=3)+
		geom_line(data=y_line_tau10bis, aes(x=x, y=y), linetype="solid", color = "#CC79A7", size=0.5)+
		geom_line(data=y_line1_tau10bis, aes(x=x, y=y), linetype="solid", color = "#E69F00", size=0.5)+
		geom_line(data=pred_lines_tau10bis, aes(x=distance, y=pdf),
				col="grey40", size=0.5, linetype="dashed")+
		scale_size_manual(name=NULL, values=c(1, rep(0.8,4), 3), guide="none")+
		scale_color_manual(name=NULL, labels=c("Great-circle KDE", "Historical (1861-1923)", "Historical (1991-2005)",
						"RCP8.5 (2006-2070)", "RCP8.5 (2071-2100)","Satellite (1982-2018)"),
				values = c("black","#00AFBB","#0072B2","#999999","#D55E00","#CC6666"),
				guide = guide_legend(override.aes = list(linetype = c(rep("solid", 5), "blank"), shape = c(rep(NA,5), 1))))+
		scale_x_log10(name="Distance (km)", limits=c(1e2, 1e5), breaks = trans_breaks("log10", function(x) 10^x),
				labels = trans_format("log10", math_format(10^.x)))+
		scale_y_log10(limits=c(10^-6, 10^-3), name = "Probability density", breaks = trans_breaks("log10", function(x) 10^x),
				labels = trans_format("log10", math_format(10^.x)))+
		theme_classic()+
		theme(legend.key.size = unit(0.15, "cm"),
				legend.position=c(0.8,0.9),
				legend.spacing = unit(0.7, "cm"))

windows(height=4, width=6)
pp_dist_tau10bis

#ggsave("~/MHWs/data/event_sync_tau10bis.pdf", useDingbats=FALSE, width = 6, height = 4)


#### ---- Supplementary figure with alternative periods for ES (tau30) ---- #####

es_dist_dat_tau30bis <- es_plotdat_tau30bis

dist.pos_tau30bis <- c(seq(1, 200, by=10), seq(201, 315, by=20), seq(316, 380, by=35),
		seq(381,460, by=40), seq(461, 800, by=50), seq(801,1100, by=70), seq(1101,2000, by=120))
dist_dat_tau30bis <- es_dist_dat_tau30bis[dist.pos_tau30bis,] %>% gather(reference, pdf, c(2:7))
dist_dat_tau30bis$line_type <- rep("s", nrow(dist_dat_tau30bis))

# to fit a power law: select a range between 10^2.6021 and 10^3.5 corresponding to 400-3162 km
power_dat_tau30bis <- es_dist_dat_tau30bis %>% filter(distance>=500&distance<=3600) %>%
		select(distance, pdf=gcdist_declust30_sat_1982_2018,
				pdf1=gcdist_declust30_1861_1923,
				pdf2=gcdist_declust30_1991_2005,
				gcdist_all)
power_fit_tau30bis <- lm(log10(pdf) ~ log10(distance), data=power_dat_tau30bis)
summary(power_fit_tau30bis)
x_range_tau30bis <- seq(500, 10^4.29, length.out=2000)
pred_lines_tau30bis <- data.frame(distance=x_range_tau30bis, pdf=10^(predict(power_fit_tau30bis, data.frame(distance=x_range_tau30bis))))

power_fit1_tau30bis <- lm(log10(pdf1) ~ log10(distance), data=power_dat_tau30bis)
summary(power_fit1_tau30bis)

power_fit2_tau30bis <- lm(log10(pdf2) ~ log10(distance), data=power_dat_tau30bis)
summary(power_fit2_tau30bis)

gc_diff_tau30bis <- es_dist_dat_tau30bis %>% filter(distance>10^3.4, distance<10^3.6) %>%
		mutate(diff=log10(gcdist_declust30_sat_1982_2018)-log10(gcdist_all)) %>%
		select(distance, diff, gcdist_all)
min.sel_tau30bis <- which(gc_diff_tau30bis$diff<0)[1]
inters_dist_tau30bis <- gc_diff_tau30bis[(min.sel_tau30bis)-1,"distance"]
inters_y_tau30bis <- es_dist_dat_tau30bis[which(es_dist_dat_tau30bis[,"distance"]>inters_dist_tau30bis)[1]-1, 2]
y_line_tau30bis <- data.frame(x=rep(inters_dist_tau30bis, 20), y=seq(0, inters_y_tau30bis, length.out=20))

gc_diff1_tau30bis <- es_dist_dat_tau30bis %>% filter(distance>10^3.6, distance<10^4) %>%
		mutate(diff=log10(gcdist_declust30_2071_2100)-log10(gcdist_all)) %>%
		select(distance, diff, gcdist_all)
min.sel1_tau30bis <- which(gc_diff1_tau30bis$diff<0)[1]
inters_dist1_tau30bis <- gc_diff1_tau30bis[(min.sel1_tau30bis)-1,"distance"]
inters_y1_tau30bis <- es_dist_dat_tau30bis[which(es_dist_dat_tau30bis[,"distance"]>inters_dist1_tau30bis)[1]-1, 2]
y_line1_tau30bis <- data.frame(x=rep(inters_dist1_tau30bis, 20), y=seq(0, inters_y1_tau30bis, length.out=20))

line_dat_tau30bis <- dist_dat_tau30bis %>% filter(reference!="gcdist_declust30_sat_1982_2018")
point_dat_tau30bis <- dist_dat_tau30bis %>% filter(reference=="gcdist_declust30_sat_1982_2018")
point_dat_tau30bis$aest <- "s"

p_dist_tau30bis <- ggplot(data=line_dat_tau30bis, aes(x=distance, y=pdf))
pp_dist_tau30bis <- p_dist_tau30bis +
		geom_line(aes(x=distance, y=pdf, col=reference, size=reference))+
		geom_point(data=point_dat_tau30bis, aes(x=distance, y=pdf, col=aest), shape=1, size=3)+
		geom_line(data=y_line_tau30bis, aes(x=x, y=y), linetype="solid", color = "#CC79A7", size=0.5)+
		geom_line(data=y_line1_tau30bis, aes(x=x, y=y), linetype="solid", color = "#E69F00", size=0.5)+
		geom_line(data=pred_lines_tau30bis, aes(x=distance, y=pdf),
				col="grey40", size=0.5, linetype="dashed")+
		scale_size_manual(name=NULL, values=c(1, rep(0.8,4), 3), guide="none")+
		scale_color_manual(name=NULL, labels=c("Great-circle KDE", "Historical (1861-1923)", "Historical (1991-2005)",
						"RCP8.5 (2006-2070)", "RCP8.5 (2071-2100)","Satellite (1982-2018)"),
				values = c("black","#00AFBB","#0072B2","#999999","#D55E00","#CC6666"),
				guide = guide_legend(override.aes = list(linetype = c(rep("solid", 5), "blank"), shape = c(rep(NA,5), 1))))+
		scale_x_log10(name="Distance (km)", limits=c(1e2, 1e5), breaks = trans_breaks("log10", function(x) 10^x),
				labels = trans_format("log10", math_format(10^.x)))+
		scale_y_log10(limits=c(10^-6, 10^-3), name = "Probability density", breaks = trans_breaks("log10", function(x) 10^x),
				labels = trans_format("log10", math_format(10^.x)))+
		theme_classic()+
		theme(legend.key.size = unit(0.15, "cm"),
				legend.position=c(0.8,0.9),
				legend.spacing = unit(0.7, "cm"))

windows(height=4, width=6)
pp_dist_tau30bis

#ggsave("~/MHWs/data/event_sync_tau30bis.pdf", useDingbats=FALSE, width = 6, height = 4)



