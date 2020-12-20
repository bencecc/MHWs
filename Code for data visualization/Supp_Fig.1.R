#### ---- Reproduce Supplementary Figure 1 ------------------------------------------ ####
#### ---- Plot maps and timeseries from network nodes. ------------------------------ ####
#### ---- Other nodes can be selected from the Shiny App at ------------------------- ####
#### ---- http://calcoloecologia.biologia.unipi.it:3838/MHW_App/   ------------------ ####
#### ---- Start selecting data for observed (global) remotely-sensed ---------------- ####
#### ---- MHWs, or historical/RCP2.6/RCP8.5 scenarios then proceed ------------------ ####
#### ---- from line 210. ------------------------------------------------------------ ####

#### ---- NOTE: here, .nc rasters are provided to color maps with the --------------- ####
#### ---- MHW "duration" covariate. Other covariates can be generated --------------- #### 
#### ---- using function net_to_map, which saves a .nc file with eight -------------- ####
#### ---- covariates: c("duration", "num_events", "intensity_cumulative", ----------- ####
#### ---- "intensity_mean", "intensity_max", "intensity_var", "rate_onset", --------- ####
#### ---- "rate_decline"). Then specify the desired variable in var ----------------- ####
#### ---- (e.g., as in line 59) for each analysis of observed or simulated MHWs. ---- ####
#### -------------------------------------------------------------------------------- ####
require(TDAmapper)
require(igraph)
#require(stars)
#require(terra)
require(sf)
require(raster)
require(rgdal)
require(ggplot2)
require(ncdf4)
require(RColorBrewer)
require(dplyr)

#load world map data as data table
load('~/MHWs/data/world_dt_robin.RData')
# source function to color nodes
source('~/MHWs/code/color_nodes_time.R')
source('~/MHWs/code/net_to_time.R')

# ----- Observed (global) --------- #
# --------------------------------- #
load('~/MHWs/data/tda_glob_res.RData')
# Load covariate data
load('~/MHWs/data/hw_glob_df.RData')
load('~/MHWs/data/tda_glob_covar.RData') #res 22, gain 45
#set .nc file to relate network nodes to map
focal_nc_obs <- '~/MHWs/data/glob_map_duration.nc'
# Generate network plot
nodes_obs <- mapperVertices(tda_glob_res, 1:13514)
edges_obs <- mapperEdges(tda_glob_res)
node_color_covar_obs <- color_nodes_time(nodes_dat=nodes_obs,
		covar_dat=tda_glob_covar, covar_name='duration')

set.seed(1567527626)
glob_graph <- graph.adjacency(tda_glob_res$adjacency, mode="undirected")
V(glob_graph)$color <- node_color_covar_obs
V(glob_graph)$node_id <- 1:nrow(nodes_obs)
my_layout_obs <- as.data.frame(layout_nicely(glob_graph, dim = 2))
my_layout_obs$node_id <- V(glob_graph)$node_id
my_layout_obs$color <- V(glob_graph)$color

sel_nodes <- c(517,524,619,620,627,726,727,734,840)
rast_dat_path <- focal_nc_obs
var <- "duration"
world_map <- world_dt_robin

# ---- Data to plot time series for selected nodes
load('~/MHWs/data/glob_years_length.RData')
# generate times for x lab
x_time_labs <- data.frame(year=glob_years_length$year,
		day=cumsum(glob_years_length$no_days), no_days=glob_years_length$no_days)
# data for plotting timeseries
pl_time_df <- data.frame(time=1:nrow(tda_glob_covar),tda_glob_covar)
nodes_sub <- nodes_obs[sel_nodes,]

# ---- Historical ----- #
# --------------------- #
sel_nodes <- c(213,214,216,218,307,308,310,311,313,396,397,398,399,400,401,494,495,496,497,
		498,499,500,622,623,624,625,626,753,755,757)
load('~/MHWs/data/tda_hist_res.RData')
# Load covariate data
load('~/MHWs/data/hw_hist_df.RData')
load('~/MHWs/data/tda_hist_covar.RData')
#set .nc file to relate network nodes to map
focal_nc_hist <- '~/MHWs/data/hist_map_duration.nc'
# Generate network plot
nodes_hist <- mapperVertices(tda_hist_res, 1:52960)
edges_hist <- mapperEdges(tda_hist_res)
node_color_covar_hist <- color_nodes_time(nodes_dat=nodes_hist,
		covar_dat=tda_hist_covar, covar_name='duration')
set.seed(1567527626)
hist_graph <- graph.adjacency(tda_hist_res$adjacency, mode="undirected")
V(hist_graph)$color <- node_color_covar_hist
V(hist_graph)$node_id <- 1:nrow(nodes_hist)

my_layout_hist <- as.data.frame(layout_nicely(hist_graph, dim = 2))
my_layout_hist$node_id <- V(hist_graph)$node_id
my_layout_hist$color <- V(hist_graph)$color

rast_dat_path <- focal_nc_hist
var <- "duration"
world_map <- world_dt_robin

# ---- Data to plot time series for selected nodes
load('~/MHWs/data/hist_years_length.RData')
# generate times for x lab
x_time_labs <- data.frame(year=hist_years_length$year,
		day=cumsum(hist_years_length$no_days), no_days=hist_years_length$no_days)
# data for plotting timeseries
pl_time_df <- data.frame(time=1:nrow(tda_hist_covar),tda_hist_covar)
nodes_sub <- nodes_hist[sel_nodes,]

# ---- RCP8.5 -------- #
# -------------------- #
load('~/MHWs/data/tda_rcp85_res.RData')
load('~/MHWs/data/hw_rcp85_df.RData')
load('~/MHWs/data/tda_rcp85_covar.RData') #res 22, gain 45

#set .nc file to relate network nodes to map
focal_nc_rcp85 <- '~/MHWs/data/rcp85_map_duration.nc'
# Generate network layout
nodes_rcp85 <- mapperVertices(tda_rcp85_res, 1:34675)
edges_rcp85 <- mapperEdges(tda_rcp85_res)
node_color_covar_rcp85 <- color_nodes_time(nodes_dat=nodes_rcp85,
		covar_dat=tda_rcp85_covar, covar_name='duration')
set.seed(1567527626)
rcp85_graph <- graph.adjacency(tda_rcp85_res$adjacency, mode="undirected")
V(rcp85_graph)$color <- node_color_covar_rcp85
V(rcp85_graph)$node_id <- 1:nrow(nodes_rcp85)

my_layout_rcp85 <- as.data.frame(layout_nicely(rcp85_graph, dim = 2))
my_layout_rcp85$node_id <- V(rcp85_graph)$node_id
my_layout_rcp85$color <- V(rcp85_graph)$color

sel_nodes <- c(877, 878, 883, 884, 888, 894, 969, 970, 971, 972, 973, 974, 976, 977, 978,
		979, 980, 983, 984, 985, 988, 989, 1058, 1059, 1060, 1061, 1062, 1063, 1065, 1066,
		1067, 1068, 1069, 1071, 1072, 1074, 1076, 1078, 1079, 1080, 1082, 1084, 1146, 1148, 1150,
		1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1161, 1162, 1165, 1166, 1167, 1168,
		1170, 1171, 1173, 1174, 1175, 1180, 1229, 1232, 1233, 1234, 1235, 1236, 1237, 1239, 1240,
		1241, 1242, 1243, 1244, 1245, 1246, 1247, 1248, 1249, 1250, 1252, 1253, 1254, 1256, 1257,
		1258, 1259, 1260, 1261, 1262, 1263, 1264, 1265, 1267, 1268, 1269, 1270, 1271, 1273, 1274,
		1275, 1276, 1277, 1278, 1279, 1280, 1281, 1282, 1283, 1284, 1285, 1286, 1288, 1289, 1290,
		1293, 1295, 1299, 1347, 1348, 1349, 1350, 1351, 1352, 1353, 1354, 1356, 1357, 1358, 1359,
		1361, 1363, 1365, 1366, 1367, 1368, 1369, 1370, 1371, 1372, 1373, 1374, 1375, 1376, 1377,
		1378, 1379, 1380, 1381, 1382, 1384, 1385, 1386, 1387, 1390, 1392, 1393, 1396, 1397, 1400,
		1403, 1445, 1446, 1447, 1448, 1449, 1450, 1451, 1453, 1454, 1456, 1457, 1458, 1459, 1460,
		1461, 1462, 1463, 1464, 1465, 1466, 1467, 1468, 1469, 1470, 1471, 1472, 1474, 1475, 1476,
		1478, 1479, 1481, 1484, 1490, 1493, 1496, 1497, 1501, 1504, 1549, 1550, 1552, 1553, 1554,
		1555, 1556, 1557, 1558, 1559, 1560, 1561, 1562, 1563, 1564, 1565, 1566, 1567, 1568, 1569,
		1570, 1571, 1572, 1573, 1575, 1576, 1577, 1580, 1581, 1582, 1583, 1584, 1585, 1586, 1591,
		1594, 1651, 1652, 1653, 1654, 1655, 1656, 1657, 1658, 1659, 1660, 1661, 1662, 1663, 1664,
		1665, 1667, 1668, 1672, 1674, 1676, 1725, 1726, 1728, 1729, 1730, 1731, 1733, 1734, 1735,
		1736, 1738, 1739, 1740, 1741, 1742, 1743, 1745, 1746, 1747, 1748, 1750, 1752, 1793, 1797,
		1798, 1802)

rast_dat_path <- focal_nc_rcp85
var <- "duration"
world_map <- world_dt_robin

# ---- Data to plot time series for selected nodes
load('~/MHWs/data/rcp85_years_length.RData')
# generate times for x lab
x_time_labs <- data.frame(year=rcp85_years_length$year,
		day=cumsum(rcp85_years_length$no_days), no_days=rcp85_years_length$no_days)
# data for plotting timeseries
pl_time_df <- data.frame(time=1:nrow(tda_rcp85_covar),tda_rcp85_covar)
nodes_sub <- nodes_rcp85[sel_nodes,]

# --------- RCP2.6 --------------#
# -------------------------------#
load('~/MHWs/data/tda_rcp26_res.RData')
# Load covariate data
load('~/MHWs/data/hw_rcp26_df.RData')
load('~/MHWs/data/tda_rcp26_covar.RData') #res 22, gain 45

#set .nc file to relate network nodes to map
focal_nc_rcp26 <- '~/MHWs/data/rcp26_map_cov.nc'
# Generate network plot
nodes_rcp26 <- mapperVertices(tda_rcp26_res, 1:34675)
edges_rcp26 <- mapperEdges(tda_rcp26_res)
node_color_covar_rcp26 <- color_nodes_time(nodes_dat=nodes_rcp26,
		covar_dat=tda_rcp26_covar, covar_name='intensity_cumulative')
set.seed(4321)
rcp26_graph <- graph.adjacency(tda_rcp26_res$adjacency, mode="undirected")
V(rcp26_graph)$color <- node_color_covar_rcp26
V(rcp26_graph)$node_id <- 1:nrow(nodes_rcp26)

my_layout_rcp26 <- as.data.frame(layout_nicely(rcp26_graph, dim = 2))
my_layout_rcp26$node_id <- V(rcp26_graph)$node_id
my_layout_rcp26$color <- V(rcp26_graph)$color

sel_nodes <- c(742,750,758,844,846,848)

rast_dat_path <- focal_nc_rcp26
var <- "duration"
world_map <- world_dt_robin

# ---- Data to plot time series for selected nodes
load('~/MHWs/data/rcp26_years_length.RData')
# generate times for x lab
x_time_labs <- data.frame(year=rcp26_years_length$year,
		day=cumsum(rcp26_years_length$no_days), no_days=rcp26_years_length$no_days)
# data for plotting timeseries
pl_time_df <- data.frame(time=1:nrow(tda_rcp26_covar),tda_rcp26_covar)
nodes_sub <- nodes_rcp26[sel_nodes,]


# ============== END DEFINING VARIABLES ===================================== #


# ------- PLOTS -------------------------------------------------------------- #

# ---- USE raster for initial data manipulations ----- #

r_tmp <- brick(rast_dat_path, varname=var)
r1 <- stack(r_tmp, layers=sel_nodes) #680,688
r1 <- flip(r1, direction=2)
if(length(sel_nodes)>1) {
	r1 <- calc(r1, fun=mean, na.rm=T)	
}

#### ---- use sf for projection and plotting (also with ggplot2) ---- ####
r1.sp <- as(r1, "SpatialPolygonsDataFrame")
r1_sf <- st_as_sf(r1.sp)
r1_sf <- r1_sf %>% st_set_crs(4326)
r1_prj <- st_transform(r1_sf, "+proj=robin")

if(names(r1_prj)[1]!="layer") {
	names(r1_prj)[1] <- "layer"
}

qn <- quantile(r1_prj$layer, c(0.99), na.rm = TRUE)	
min_val <- ifelse(min(r1_prj$layer, na.rm=T)>=0,min(r1_prj$layer, na.rm=T),0)
max_val <- max(r1_prj$layer, na.rm=T)
range_val <- range(r1_prj$layer, na.rm=T)
max_val_below_qn <- as.numeric(round(data.frame(value=r1_prj$layer) %>% filter(value<=qn) %>% summarise(max_val_bqn=max(value)), digits=1))

col_pal <- rev(colorRampPalette(brewer.pal(11, 'Spectral'))(14))

resp_var_name <- switch(var,
		"duration"="Duration (days)",
		"num_events"="Number of events",
		"intensity_cumulative"=expression(atop("Cumulative intensity","("*degree*C~"x days)")),
		"intensity_mean"=expression(paste("Mean intensity ","("*degree*C~")")),
		"intensity_max"=expression(paste("Maximum intensity ","("*degree*C~")")),
		"intensity_var"=expression(paste("Variance intensity (",degree,"C"^2,")")),
		"rate_onset"=expression(atop("Rate onset","("*degree*C~"x days"^-1*")")),
		"rate_decline"=expression(atop("Rate decline","("*degree*C~"x days"^-1*")")))

# ----- Bounding box (frame - see below to generate a polygon with sf ----- #
len.s <- 1000
frame_robin <- data.frame(
		x=c(seq(-180,180,length.out = len.s),
				rep(-180,len.s),
				seq(180,-180,length.out = len.s),
				rep(180,len.s)),
		y=c(rep(-90, len.s),
				seq(-90,90,length.out = len.s),
				rep(90, len.s),
				seq(90,-90,length.out = len.s)))


data.dum <- data.frame(path=rep(1,nrow(frame_robin)))
robin_frm_tmp <- st_as_sf(SpatialPointsDataFrame(coords = frame_robin,
				data=data.dum))
robin_frm <- robin_frm_tmp %>% st_set_crs(4326) %>% st_transform(crs="+proj=robin")
robin_df <- data.frame(st_coordinates(robin_frm))
robin_df$path <- robin_frm$path

hwplot <- ggplot()+
#		geom_sf(data=r1_prj, aes(col=layer))+
		geom_sf(data = r1_prj %>% filter(layer<=qn), aes(fill=layer, col=layer)) + 
		geom_sf(data = r1_prj %>% filter(layer>qn), fill=col_pal[14], col=col_pal[14]) +
		scale_color_gradientn(name=resp_var_name,
				colours = col_pal,
				guide = "colorbar", na.value='green',
				labels=scales::scientific) +
		scale_fill_gradientn(name=resp_var_name,
				colours = col_pal,
				guide = "colorbar", na.value='green',
				labels=scales::scientific) +
		geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
				fill="black", col="black", size = 0.5) +
		geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
				fill="gray80") +
		geom_path(data = robin_df, aes(x=X, y=Y), col="grey40",
				alpha=1, size=0.2) +
		coord_sf(datum=NA, ndiscr=100) +
		theme_void() +
		theme(legend.text=element_text(size=8),
				legend.title=element_text(size=8),
				legend.spacing.x = unit(0.1, 'cm'),
				legend.position = 'right',
				legend.key.size = unit(0.6, "cm"),
				plot.margin=unit(c(-1,0,0,0),"cm"))

windows(height=5, width=6)
plot(hwplot)

#ggsave("~/Lavori/ExtremeEvents/CMIP5/Plots/glob_map_blob.pdf", useDingbats=FALSE, width = 6, height = 3)
#ggsave("~/Lavori/ExtremeEvents/CMIP5/Plots/hist_map_2001-2003.pdf", useDingbats=FALSE, width = 6, height = 3)
#ggsave("~/Lavori/ExtremeEvents/CMIP5/Plots/rcp26_map_2099-2100.pdf", useDingbats=FALSE, width = 6, height = 3)
#ggsave("~/Lavori/ExtremeEvents/CMIP5/Plots/rcp85_map_2068-2100.pdf", useDingbats=FALSE, width = 6, height = 3)


# ---- Plot time series for selected nodes
focal_tmp <- apply(nodes_sub, 1, function(x) {
			str.tmp <- strsplit(as.character(x['Nodename']), split=c('V',':'))[[1]][2]
			substr(str.tmp, start=regexpr(':', str.tmp)[1], stop=regexpr(':', str.tmp)[1]) <- ','
			focal_times=as.vector(as.numeric(unlist(strsplit(str.tmp,","))))[-1]
		})	

#unique times
if(is.list(focal_tmp)) {
	focal_times <- data.frame(time=unlist(focal_tmp)) %>% distinct() %>% arrange(time)
}

if(!is.list(focal_tmp)) {
	focal_times <- data.frame(time=as.vector(focal_tmp)) %>% distinct() %>% arrange(time)
}

if(nrow(focal_times)==1) {
	
	t_sel <- x_time_labs %>% filter(day>=focal_times$time) %>% slice(., 1L)
	
	if(t_sel$day==focal_times$time) {
		day_sel <- t_sel$day
		year_sel <- droplevels(t_sel$year)
	}
	
	else {
		row_sel <- which(x_time_labs$day%in%t_sel$day)-1
		day_sel <- focal_times$time-x_time_labs[row_sel,'day']
		year_sel <- droplevels(x_time_labs[row_sel,'year'])
	}
	
	xlab <- paste('day ', day_sel, ' of year ', year_sel, sep="")
	
	pl_t <- ggplot(data=data.frame(pl_time_df[focal_times$time,]), aes_string(x='time', y=var))+
			geom_col(aes(fill=var), width=0.5, position='identity', show.legend=F)+
			scale_x_discrete(labels=xlab)+
			ggtitle('Only one day selected')+
			labs(x=xlab, y=resp_var_name)+
			theme_bw(base_size = 12)+
			theme(panel.grid.major=element_line(colour=NA),
					axis.text.x=element_blank(),
					axis.title.y=element_text(size=11),
					panel.grid.minor=element_line(colour=NA),
					plot.title = element_text(hjust = 0.5))
	
}

if(nrow(focal_times)>1) {
	
	t_sel <- x_time_labs %>% filter(day%in%focal_times$time)
	range_time <- range(focal_times$time)
	t_sel <- x_time_labs %>% filter(day>=range_time[1], day<=range_time[2]) #do(data.frame(sel=intersect(.$day, range_time[1]:range_time[2])))
	
	# if times are within a single year
	if(nrow(t_sel)==0){
		t_sel <- x_time_labs %>% filter(day>=range_time[1]) %>% slice(., 1L)
		row_sel <- which(x_time_labs$day%in%t_sel$day)-1
		
		if(row_sel==0) {
			row_sel=1
			tick_lab <- focal_times$time
		}
		
		else {
			tick_lab <- focal_times$time-x_time_labs[row_sel,'day']
		}
		
		year_sel <- droplevels(x_time_labs[row_sel,'year'])
		xlab <- paste('days in year ', year_sel, sep="")
		# get days to label x axis
		int_time <- focal_times$time
		angle_sel <- 0
	}
	
	else {
		
		# get years to label x axis
		tick_lab <- droplevels(t_sel$year)
		int_time <- t_sel$day
		xlab <- ""
		angle_sel <- 45
	}
	
	#reduce number of ticks if necessary
	if(length(int_time)>20) {
		
		if(length(int_time)<=50) {
			int_time <- int_time[c(TRUE,FALSE)]
			tick_lab <- tick_lab[c(TRUE,FALSE)]
		}
		
		if (length(int_time)>50&&length(int_time)<=100) {
			int_time <- int_time[c(TRUE,rep(FALSE, 5))]
			tick_lab <- tick_lab[c(TRUE,rep(FALSE, 5))]
		}
		
		if (length(int_time)>100) {
			int_time <- int_time[c(TRUE,rep(FALSE, 11))]
			tick_lab <- tick_lab[c(TRUE,rep(FALSE, 11))]
		}
		
		angle_sel <- 90
	}
	
	pl_t <- ggplot(data=pl_time_df[focal_times$time,], aes_string(x='time', y=var))+
#			geom_path(color='blue')+
			geom_path(color="#FC4E07", size=1)+
#			geom_point(fill="#4E84C4", col="#4E84C4", shape=21, size=2, show.legend=F)+
#			geom_point(col="#FC4E07", fill="#FC4E07", shape=21, size=1.5, show.legend=F)+
			scale_x_continuous(breaks=int_time, expand=c(0.05,0.05), 
					labels=tick_lab)+ 
			labs(x=xlab, y=resp_var_name)+
			theme_bw(base_size = 12)+
			theme(panel.grid.major=element_line(colour=NA),
					axis.text.x=element_text(size=10, angle=45, vjust=0.7),
					axis.text.y=element_text(size=10, angle=0),
					axis.title.x=element_text(size=16),
					axis.title.y=element_text(size=16),
					panel.grid.minor=element_line(colour=NA),
					plot.title = element_text(hjust = 0.5))
	#plot.margin=unit(c(-1,0,0,0),"cm"))
}


windows(width=6, height=3)
plot(pl_t)

#ggsave("~/Lavori/ExtremeEvents/CMIP5/Plots/timeser_glob_map_blob.pdf", useDingbats=FALSE, width = 6, height = 3)
#ggsave("~/Lavori/ExtremeEvents/CMIP5/Plots/timeser_hist_map_2001-2003.pdf", useDingbats=FALSE, width = 6, height = 3)
#ggsave("~/Lavori/ExtremeEvents/CMIP5/Plots/timeser_rcp26_map_2099-2100.pdf", useDingbats=FALSE, width = 6, height = 3)
#ggsave("~/Lavori/ExtremeEvents/CMIP5/Plots/timeser_rcp85_map_2068-2100.pdf", useDingbats=FALSE, width = 6, height = 3)









