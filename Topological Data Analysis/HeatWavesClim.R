# The function detects heatwaves from .nc files and generates data structure
# for input # in TDA analysis. Here, a climatology is obtained from historical
# data. The function is called by TDA_pipeline.R
###############################################################################
#!/bin/bash Rscript

require(ncdf4)
require(heatwaveR)
require(ncdf4.helpers)
library(PCICt)
require(raster)
require(dplyr)

require(foreach, quietly=T)
require(doMC, quietly=T)
registerDoMC(cores=20)

# set focal_var: tos for CMIP data, sst for SSTglobal data
focal_var <- 'tos'

args <- commandArgs(trailingOnly = TRUE)
dat_nc <- nc_open(args[1])
lat_first <- as.numeric(args[2])

cat(lat_first)

##### ---- Generate data structure to compute climatology from historical model data ---- ####
dat_hist <- nc_open('~/MHWs/data/hist.ens.nc') # hist_0.nc contains the 1861-2005 historical data

get_clim <- function(dat_nc, which_lat, ...) {
	
	tos_time <- nc.get.time.series(dat_nc, v = "tos",
			time.dim.name = "time")
	timestamp = as.Date(format(tos_time, "%Y-%m-%d"))
	dat <- ncvar_get(dat_nc, varid="tos", start=c(1,which_lat,1), count=c(-1,1,-1))
	list(dat,timestamp)
}

dat_clim <- get_clim(dat_nc=dat_hist, which_lat=lat_first)
rm(dat_hist)
###############################################################################################

#focal_var is tos for CMIP data, sst for SSTglobal data
tos_time <- nc.get.time.series(dat_nc, v = focal_var,
		time.dim.name = "time")
timestamp = as.Date(format(tos_time, "%Y-%m-%d"))

years_all <- data.frame(year=format(tos_time, "%Y"))
years_length <- years_all %>% group_by(year) %>% summarise(no_days=n())

focal_doy <- NULL
for (i in 1:nrow(years_length)) {
	if(years_length[i,'no_days']==365) {
		focal_doy <- c(focal_doy, c(1:59,61:366))
	}
	else {
		focal_doy <- c(focal_doy, 1:366)
	}	
}

dat <- ncvar_get(dat_nc, varid=focal_var, start=c(1,lat_first,1), count=c(-1,1,-1))
long.vec <- ncvar_get(dat_nc,"lon")
lat.vec <- ncvar_get(dat_nc,"lat")

out <- foreach(i=1:nrow(dat), .combine='rbind',
				.packages=c('ncdf4','heatwaveR'),
				.export=c('dat', 'lat_first')) %dopar% {
			
			long <- long.vec[i]
			lat <- lat.vec[lat_first]
			
			lat_id <- ifelse(lat>0, abs(lat-90.5), 90.5+abs(lat))
			
			if ((lat_id-1)*360+long+1<=360) {
				cell_id <- long+1
			}
			
			else {
				cell_id <- (lat_id-1)*360+long+1
			}
			
			exc <- which(is.na(dat[i,]))
			exc1 <- which(is.na(dat_clim[[1]][i,]))
			
			# select only time series that have less than 10% of NAs
			if (length(exc)>0.1*ncol(dat)|length(exc1)>0.1*ncol(dat_clim[[1]])) {
				res <- NULL
			}
			
			else {
				test_dat <- data.frame(t=as.Date(as.character(dat_clim[[2]])), temp=(dat_clim[[1]][i,]-273.15))
				ts <- try(ts2clm(test_dat, climatologyPeriod = c("1861-01-01", "2005-12-31")),
						silent=T)
				
				if (inherits(ts, "try-error")) {
					res <- NULL
				}
				
				else  {	
					
					#### ---- Uncomment lines below to use with rcp85 and rcp26 ---- #####
					# reference data for climatology; variables of interest are seas and thresh
					# to use with rcp model data; julian day 60 is removed because rcp8.5 and rcp2.6 models
					# do not include leap years
					ref_clim <- ts %>% distinct(doy, .keep_all=T) %>% filter(doy!=60)
					
					# to get seas and thresh to match model data it is necessary to duplicate 
					# the (julian year) ref_clim object for the number of model years (95, from 2006  to 2100)
					
					no_model_yrs <- length(focal_doy)/365	
					
					new_ts <- data.table(doy=focal_doy, t=as.Date(as.character(timestamp)),
							temp=dat[i,]-273.15, seas=rep(ref_clim$seas, no_model_yrs),
							thresh=rep(ref_clim$thresh,no_model_yrs))
########################################################################
					
#### ---- Uncomment lines below to use with SSTglobal data ---- ####
#		#find leap years
#      leap <- which(ts$doy==60)
#      #select climatolgy for the 366 days in a leap year
#      ts_366 <- ts[(leap[1]-59):(leap[1]+306),]
#      #select climatology for the 365 days in a normal year
#      ts_365 <- ts[1:365,]
#      
#      #generate seas and thresh vectors with appropriate climatology for leap and not-leap years 
#      seas <- NULL
#      thresh <- NULL
#      
#      for (j in 1:nrow(years_length)) {
#        
#        if (years_length[j,'no_days']==365) {
#          seas <- c(seas, ts_365$seas)
#          thresh <- c(thresh, ts_365$thresh)
#          }
#        
#        else {
#          seas <- c(seas, ts_366$seas)
#          thresh <- c(thresh, ts_366$thresh)
#        }
#                  
#      }
#      
#      new_ts <- data.table(doy=focal_doy, t=as.Date(as.character(timestamp)),
#					temp=dat[i,], seas=seas, thresh=thresh)
###########################################################################
					
					hw.res <- try(detect_event(new_ts), silent=T)
					
					if (inherits(hw.res, "try-error")|length(hw.res$event$event_no)==0) {
						res <- NULL
					}
					
					else  {	      
						start_date <- round(hw.res$event$index_start, 1) 
						cx <- rep(long, length(start_date))
						cy <- rep(lat, length(start_date))
						res <- data.frame(cell=cell_id, x=cx, y=cy, hw.res$event)
					}
				}
			}
			
			plyr::compact(res)
			
		}

nc_close(dat_nc)

outputName=paste('lat_', lat_first, ".RData", sep="")
outputPath=file.path('~/tmp_files', outputName)
assign(paste('lat_', lat_first, sep=""), value=out, pos=1, inherits=T)
save(list=paste('lat_', lat_first, sep=""), file=outputPath)



