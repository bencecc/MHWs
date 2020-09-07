# The function detects heatwaves from .nc files and generates data structure
# for input in TDA analysis. The function derives a climatology from tis own
# data; see # function HeatWavesClim.R to derive a climatology from a different
# dataset # (e.g., historical climatology). The function is called by
# TDA_pipeline.R
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

# set focal_var: tos (in kelvin) for CMIP data, sst for SSTglobal data
focal_var <- 'sst'

args <- commandArgs(trailingOnly = TRUE)
dat_nc <- nc_open(args[1])
lat_first <- as.numeric(args[2])
cat(lat_first)

sst_time <- nc.get.time.series(dat_nc, v = focal_var,
		time.dim.name = "time")
timestamp = as.Date(format(sst_time, "%Y-%m-%d"))

years_all <- data.frame(year=format(sst_time, "%Y"))
years_length <- years_all %>% group_by(year) %>% summarise(no_days=n())

doy <- NULL
for (i in 1:nrow(years_length)) {
	if(years_length[i,'no_days']==365) {
		doy <- c(doy, c(1:59,61:366))
	}
	else {
		doy <- c(doy, 1:366)
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
			
			if (length(exc)>0.1*ncol(dat)) {
				res <- NULL
			}
			
			else {
				
				test_dat <- data.frame(t=as.Date(as.character(timestamp)), temp=dat[i,]) #-273.15 if variable tos in kelvin
				ts <- try(ts2clm(test_dat, climatologyPeriod = c("1982-01-01", "2018-12-31")), 
						#climatology is c("1861-01-01", "2005-12-31") for historical data "1982-01-01", "2018-12-31" for contemporary data
						silent=T)
				
				if (inherits(ts, "try-error")) {
					res <- NULL
				}
				
				else  {	
					
					hw.res <- try(detect_event(ts), silent=T)
					
					if (inherits(hw.res, "try-error")) {
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
