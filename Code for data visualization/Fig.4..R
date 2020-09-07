#### ---- Distance distribution of significantly synchronized connection - Fig. 4 ---- ####
#### --------------------------------------------------------------------------------- ####

require(ggplot2)
require(scales)
require(dplyr)
require(tidyr)

load("~/MHWs/data/es_plotdat_tau10.RData")

#### ---- Plot Fig 4 of main manuscript (tau = 10) ---- ####
es_dist_dat_tau10 <- es_plotdat_tau10

dist.pos_tau10 <- c(seq(1, 200, by=20), seq(201, 315, by=40), seq(316, 380, by=80),
		seq(381,560, by=100), seq(561,1100, by=80), seq(1101,2000, by=120))
dist_dat_tau10 <- es_dist_dat_tau10[dist.pos_tau10,] %>% gather(reference, pdf, c(2:7))
dist_dat_tau10$line_type <- rep("s", nrow(dist_dat_tau10))

# to fit a power law: select a range between 10^2.6021 and 10^3.5 corresponding to 400-3162 km
power_dat_tau10 <- es_dist_dat_tau10 %>% filter(distance>=500&distance<=3600) %>%
		select(distance, pdf=gcdist_declust_sat_1982_2018,
				pdf1=gcdist_declust_1891_1900,
				pdf2=gcdist_declust_1991_2000,
				gcdist_all)
power_fit_tau10 <- lm(log10(pdf) ~ log10(distance), data=power_dat_tau10)
summary(power_fit_tau10)
x_range_tau10 <- seq(500, 10^4.29, length.out=2000)
pred_lines_tau10 <- data.frame(distance=x_range_tau10, pdf=10^(predict(power_fit_tau10, data.frame(distance=x_range_tau10))))

power_fit1_tau10 <- lm(log10(pdf1) ~ log10(distance), data=power_dat_tau10)
summary(power_fit1_tau10)

power_fit2_tau10 <- lm(log10(pdf2) ~ log10(distance), data=power_dat_tau10)
summary(power_fit2_tau10)

gc_diff_tau10 <- es_dist_dat_tau10 %>% filter(distance>10^3.4, distance<10^3.6) %>%
		mutate(diff=log10(gcdist_declust_sat_1982_2018)-log10(gcdist_all)) %>%
		select(distance, diff, gcdist_all)
min.sel_tau10 <- which(gc_diff_tau10$diff<0)[1]
inters_dist_tau10 <- gc_diff_tau10[(min.sel_tau10)-1,"distance"]
inters_y_tau10 <- es_dist_dat_tau10[which(es_dist_dat_tau10[,"distance"]>inters_dist_tau10)[1]-1, 2]
y_line_tau10 <- data.frame(x=rep(inters_dist_tau10, 20), y=seq(0, inters_y_tau10, length.out=20))

gc_diff1_tau10 <- es_dist_dat_tau10 %>% filter(distance>10^3.6, distance<10^3.8) %>%
		mutate(diff=log10(gcdist_declust_2081_2090)-log10(gcdist_all)) %>%
		select(distance, diff, gcdist_all)
min.sel1_tau10 <- which(gc_diff1_tau10$diff<0)[1]
inters_dist1_tau10 <- gc_diff1_tau10[(min.sel1_tau10)-1,"distance"]
inters_y1_tau10 <- es_dist_dat_tau10[which(es_dist_dat_tau10[,"distance"]>inters_dist1_tau10)[1]-1, 2]
y_line1_tau10 <- data.frame(x=rep(inters_dist1_tau10, 20), y=seq(0, inters_y1_tau10, length.out=20))

line_dat_tau10 <- dist_dat_tau10 %>% filter(reference!="gcdist_declust_sat_1982_2018")
point_dat_tau10 <- dist_dat_tau10 %>% filter(reference=="gcdist_declust_sat_1982_2018")
point_dat_tau10$aest <- "s"

p_dist_tau10 <- ggplot(data=line_dat_tau10, aes(x=distance, y=pdf))
pp_dist_tau10 <- p_dist_tau10 +
		geom_line(aes(x=distance, y=pdf, col=reference, size=reference))+
		geom_point(data=point_dat_tau10, aes(x=distance, y=pdf, col=aest), shape=1, size=3)+
		geom_line(data=y_line_tau10, aes(x=x, y=y), linetype="solid", color = "#CC79A7", size=0.5)+
		geom_line(data=y_line1_tau10, aes(x=x, y=y), linetype="solid", color = "#E69F00", size=0.5)+
		geom_line(data=pred_lines_tau10, aes(x=distance, y=pdf),
				col="grey40", size=0.5, linetype="dashed")+
		scale_size_manual(name=NULL, values=c(1, rep(0.8,4), 3), guide="none")+
		scale_color_manual(name=NULL, labels=c("Great-circle KDE", "Historical (1891-1900)", "Historical (1991-2000)",
						"RCP8.5 (2041-2050)", "RCP8.5 (2081-2090)","Satellite (1982-2018)"),
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
pp_dist_tau10

#ggsave("~/MHWs/data/event_sync_tau10.pdf", useDingbats=FALSE, width = 6, height = 4)


