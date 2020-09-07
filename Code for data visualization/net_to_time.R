#### ---- Function to display time series of node ------- ####
#### ---- attributes from selected nodes ina network ---- ####
#### ---------------------------------------------------- ####

net_to_time <- function(dat, focal_nodes=NULL,
		plot.var, time_sel, plot=F, ...) {

	ylab <- switch(plot.var,
			"intensity_cumulative"=expression(atop("intensity_cumulative","("*degree*C~"x days)")),
			"duration"="duration (days)",
			"intensity_max"="intensity_max (°C)",
			"intensity_mean"="intensity_mean (°C)",
			"intensity_max"="intensity_max (°C)",
			"intensity_var"=expression(paste("intensity_var (",~degree,"C"^2," )")),
			"rate_onset"=expression(atop("rate_onset","("*degree*C~"x days"^-1*")")),
			"rate_decline"=expression(atop("rate_decline","("*degree*C~"x days"^-1*")")),
	)
		
	focal_tmp <- apply(focal_nodes, 1, function(x) {
				str.tmp <- strsplit(as.character(x['Nodename']), split=c('V',':'))[[1]][2]
				substr(str.tmp, start=regexpr(':', str.tmp)[1], stop=regexpr(':', str.tmp)[1]) <- ','
				focal_times=as.vector(as.numeric(unlist(strsplit(str.tmp,","))))[-1]
			})	
	
	#unique times
	if(is.list(focal_tmp)) {
		focal_times <- data.frame(time=unlist(focal_tmp)) %>% distinct() %>% arrange(time)
	}
	
	else {
		focal_times <- data.frame(time=as.vector(focal_tmp)) %>% distinct() %>% arrange(time)
	}
	
	if(nrow(focal_times)==1) {
		
		t_sel <- time_sel %>% filter(day>=focal_times$time) %>% slice(., 1L)
		
		if(t_sel$day==focal_times$time) {
			day_sel <- t_sel$day
			year_sel <- droplevels(t_sel$year)
		}
		
		else {
			row_sel <- which(time_sel$day%in%t_sel$day)-1
			day_sel <- focal_times$time-time_sel[row_sel,'day']
			year_sel <- droplevels(time_sel[row_sel,'year'])
		}
		
		xlab <- paste('day ', day_sel, ' of year ', year_sel, sep="")
		
		pl_t <- ggplot(data=data.frame(dat[focal_times$time,]), aes_string(x='time', y=plot.var))+
				geom_col(aes(fill=plot.var), width=0.5, position='identity', show.legend=F)+
				scale_x_discrete(labels=xlab)+
				ggtitle('Only one day selected')+
				labs(x=xlab, y=ylab)+
				theme_bw(base_size = 12)+
				theme(panel.grid.major=element_line(colour=NA),
						axis.text.x=element_blank(),
						axis.title.y=element_text(size=11),
						panel.grid.minor=element_line(colour=NA),
						plot.title = element_text(hjust = 0.5))
		
	}
	
	if(nrow(focal_times)>1) {
		
		t_sel <- time_sel %>% filter(day%in%focal_times$time)
		range_time <- range(focal_times$time)
		t_sel <- time_sel %>% filter(day>=range_time[1], day<=range_time[2]) #do(data.frame(sel=intersect(.$day, range_time[1]:range_time[2])))
		
		# if times are within a single year
		if(nrow(t_sel)==0){
			t_sel <- time_sel %>% filter(day>=range_time[1]) %>% slice(., 1L)
			row_sel <- which(time_sel$day%in%t_sel$day)-1
			
			if(row_sel==0) {
				row_sel=1
				tick_lab <- focal_times$time
			}
			
			else {
				tick_lab <- focal_times$time-time_sel[row_sel,'day']
			}
			
			year_sel <- droplevels(time_sel[row_sel,'year'])
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
		
		pl_t <- ggplot(data=dat[focal_times$time,], aes_string(x='time', y=plot.var))+
				geom_point(aes(col=plot.var), show.legend=F)+
				geom_path(color='blue')+
				scale_x_continuous(breaks=int_time, expand=c(0.05,0.05), 
						labels=tick_lab)+ 
				ggtitle('Timeseries from selected nodes')+
				labs(x=xlab, y=ylab)+
				theme_bw(base_size = 12)+
				theme(panel.grid.major=element_line(colour=NA),
						axis.text.x=element_text(size=9, angle=angle_sel, hjust=1),
						axis.title.y=element_text(size=11),
						panel.grid.minor=element_line(colour=NA),
						plot.title = element_text(hjust = 0.5))
						
	}
	
	if(plot){
		plot(pl_t)
	}
	
	else {
		pl_t
	}
}

