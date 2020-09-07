#### ---- Function to generate covariates for quick annotation of maps from selected ----------- ####
#### ---- network nodes. It aggregates MHW statistics by node instead of time, ----------------- ####
#### ---- although the third variable is called Time. Use in parallel or with multithreads; ---- ####
#### ---- file_init is a character giving the prefix for saved .nc files; hw_res is a df ------- ####
#### ---- with info on MHW such as hw_hist_df; net_nodes are the nodes to be mapped as --------- ####
#### ---- a function of one of the coloring variables. ----------------------------------------- ####
#### ------------------------------------------------------------------------------------------- ####

require(TDAmapper)
require(dplyr)
require(ncdf4)
require(ncdf4.helpers)
require(easyNCDF)
require(raster)

require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=70)

net_to_map <- function(file_init, hw_res, net_nodes, ...) {
	
	if(!file.exists('~/tmp_map_cov'))
		dir.create('~/tmp_map_cov')
	
	setwd('~/tmp_map_cov')
	
	tmp_r <- raster(nrows=180, ncols=360)
	xmin(tmp_r) <- -180
	xmax(tmp_r) <- 180
	ymin(tmp_r) <- -90
	ymax(tmp_r) <- 90
	
	# get days within nodes
	focal_times <- apply(net_nodes, 1, function(x) {
				str.tmp <- strsplit(as.character(x['Nodename']), split=c('V',':'))[[1]][2]
				substr(str.tmp, start=regexpr(':', str.tmp)[1], stop=regexpr(':', str.tmp)[1]) <- ','
				focal_times=as.vector(as.numeric(unlist(strsplit(str.tmp,","))))[-1]
			})	
	
	# generate ncdf file with variables reflecting aggregated (mean, sum, total)
	# MHW properties within each node
	out_list <- foreach(i=1:length(focal_times)) %dopar% {
		
		cat('Doing iter ', i, ' of ', length(focal_times), '\n', sep = '')
		
		t <- ncdim_def( "Time", "days since 2006-01-01", i, unlim=TRUE) # Time is actually node number
		x <- ncdim_def( "lon", "degreesE", -179.5:179.5)
		y <- ncdim_def( "lat", "degreesN", as.double(-89.5:89.5))
		
		n1 <- ncvar_def("duration", "days",
				list(x,y,t),-9999, prec = "double")
		n2 <- ncvar_def("num_events", "number_of_mhws",
				list(x,y,t),-9999, prec = "double")
		n3 <- ncvar_def("intensity_cumulative", "deg_C x days",
				list(x,y,t),-9999, prec = "double")
		n4 <- ncvar_def("intensity_mean", "deg_C",
				list(x,y,t),-9999, prec = "double")
		n5 <- ncvar_def("intensity_max", "deg_C",
				list(x,y,t),-9999, prec = "double")
		n6 <- ncvar_def("intensity_var", "deg_C^2",
				list(x,y,t),-9999, prec = "double")
		n7 <- ncvar_def("rate_onset", "deg_C x days^-1",
				list(x,y,t),-9999, prec = "double")
		n8 <- ncvar_def("rate_decline", "deg_C x days^-1",
				list(x,y,t),-9999, prec = "double")
		
		
		# each element of this list is a tibble with MHW data for each
		# of the days included in the node
		out_tmp <- list()
		
		for(j in 1:length(focal_times[[i]])) {
			
			cross_time <- focal_times[[i]][j]
			
			node_geo <- hw_res %>%
					do(data.frame(sel=data.table::between(cross_time, .$index_start, .$index_end))) %>%
					do(data.frame(sel=which(.$sel==T))) %>%
					do(data.frame(hw_res[.$sel,])) 
			
			node_geo$x <- node_geo$x-179.5
			shift <- 20
			node_geo$x <- node_geo$x+shift
			node_geo$x <- ifelse(node_geo$x > 179.5, node_geo$x-360, node_geo$x)
			
			out_tmp[[j]] <- node_geo 
		}
		
		# MHW statistics grouped by lon and lat and aggregated across days in the node
		tmp_list <- bind_rows(out_tmp) %>% group_by(x,y) %>% distinct(event_no, .keep_all=T) %>% #distinct prevents counting the same event more than once in a node
				summarise(duration=mean(duration, na.rm=T),
						num_events=n(), 
						intensity_cumulative=sum(intensity_cumulative, na.rm=T), #sum across days in a node
						intensity_mean=mean(intensity_mean, na.rm=T),
						intensity_max=max(intensity_max, na.rm=T),
						intensity_var=mean(intensity_var, na.rm=T),
						rate_onset=mean(rate_onset, na.rm=T),
						rate_decline=mean(rate_decline, na.rm=T)) %>% arrange(desc(y), x)
	
		attr(tmp_list, 'indices') <- attr(tmp_list, 'vars') <- attr(tmp_list, 'drop') <- 
				attr(tmp_list, 'group_sizes') <- attr(tmp_list, 'biggest_group_size') <- attr(tmp_list, 'labels') <- NULL
		
		file_save <- sprintf(paste('%04d.nc', sep=""),i)
		ncout <- nc_create(file_save, list(n1,n2,n3,n4,n5,n6,n7,n8), force_v4=T)
		var_names <- NcReadVarNames(file_save)[1:8]
		
		r1 <- rasterize(as.matrix(tmp_list[,c('x','y')]), tmp_r, as.matrix(tmp_list[,'duration']), fun=mean, na.rm=T)
		ncvar_put(ncout, n1, values(r1), start=c(1,1,1), count=c(-1,-1,1))
		r2 <- rasterize(as.matrix(tmp_list[,c('x','y')]), tmp_r, as.matrix(tmp_list[,'num_events']), fun=mean, na.rm=T)
		ncvar_put(ncout, n2, values(r2), start=c(1,1,1), count=c(-1,-1,1))
		r3 <- rasterize(as.matrix(tmp_list[,c('x','y')]), tmp_r, as.matrix(tmp_list[,'intensity_cumulative']), fun=mean, na.rm=T)
		ncvar_put(ncout, n3, values(r3), start=c(1,1,1), count=c(-1,-1,1))
		r4 <- rasterize(as.matrix(tmp_list[,c('x','y')]), tmp_r, as.matrix(tmp_list[,'intensity_mean']), fun=mean, na.rm=T)
		ncvar_put(ncout, n4, values(r4), start=c(1,1,1), count=c(-1,-1,1))
		r5 <- rasterize(as.matrix(tmp_list[,c('x','y')]), tmp_r, as.matrix(tmp_list[,'intensity_max']), fun=mean, na.rm=T)
		ncvar_put(ncout, n5, values(r5), start=c(1,1,1), count=c(-1,-1,1))
		r6 <- rasterize(as.matrix(tmp_list[,c('x','y')]), tmp_r, as.matrix(tmp_list[,'intensity_var']), fun=mean, na.rm=T)
		ncvar_put(ncout, n6, values(r6), start=c(1,1,1), count=c(-1,-1,1))
		r7 <- rasterize(as.matrix(tmp_list[,c('x','y')]), tmp_r, as.matrix(tmp_list[,'rate_onset']), fun=mean, na.rm=T)
		ncvar_put(ncout, n7, values(r7), start=c(1,1,1), count=c(-1,-1,1))
		r8 <- rasterize(as.matrix(tmp_list[,c('x','y')]), tmp_r, as.matrix(tmp_list[,'rate_decline']), fun=mean, na.rm=T)
		ncvar_put(ncout, n8, values(r8), start=c(1,1,1), count=c(-1,-1,1))
		
		ncatt_put(ncout,"Time","axis","T") # Time is actually node number
		ncatt_put(ncout,"lon","axis","X")
		ncatt_put(ncout,"lat","axis","Y")
		
		nc_close(ncout)
		
	}
	
}

load('~/MHWs/data/hw_glob_df.RData')
load('~/MHWs/data/tda_glob_res.RData')
nodes <- mapperVertices(tda_glob_res, 1:13514)

glob_map_cov <- net_to_map(file_init="", hw_res=hw_glob_df,
		net_nodes=nodes)

setwd('~/tmp_map_cov') 
system("cdo cat *.nc glob_map_cov.nc")


