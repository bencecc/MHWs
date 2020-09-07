#### ---- This function aggregates covariate values based on the elements ------ ####
#### ---- in a node and allows the annotation of nodes in a shape graph; ------- ####
#### ---- it takes as arguments a mapperVertices object and is used on wide ---- ####
#### ---- (times as rows, cells as columns) data. Argument include: ------------ ####
#### ---- nodes_dat, a dataframe containing informaiton on network nodes, ------ ####
#### ---- as generated from TDA function MapperVertices; covar_dat, dataset ---- ####
#### ---- of covariates (e.g., tda_glob_covar) and covar_name (name of the ----- ####
#### ---- covariate used to annotate the network graph. ------------------------ ####
#### --------------------------------------------------------------------------- ####

color_nodes_time <- function(nodes_dat, covar_dat, covar_name, ...) { 
	
	apply(nodes_dat, 1, function(x) {
				
				#identify the rows of the time_dat contributing to the node
				str.tmp <- strsplit(as.character(x['Nodename']), split=c('V',':'))[[1]][2]
				substr(str.tmp, start=regexpr(':', str.tmp)[1], stop=regexpr(':', str.tmp)[1]) <- ','
				focal_times=as.vector(as.numeric(unlist(strsplit(str.tmp,","))))[-1]
				
				#select times on the covar dataset
				sub_dat <- covar_dat[focal_times,]
				
				if (covar_name=='num_events'|covar_name=='intensity_cumulative') {
					sum(sub_dat[, covar_name], na.rm=T)
				}
				
				else if (covar_name=='intensity_max') {
					max(sub_dat[, covar_name], na.rm=T)
				}
				
				else if (covar_name=='duration'|covar_name=='intensity_mean'|covar_name=='intensity_var'
						|covar_name=='rate_onset'|covar_name=='rate_decline') {
					mean(sub_dat[, covar_name], na.rm=T)
				}
			})
}




		
		