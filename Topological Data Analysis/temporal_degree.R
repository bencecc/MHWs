# Derivation of node degree from a tda_res object obtained as an output of a
# TDA analysis (e.g., tda_optim). The function is called by TDA_pipeline.R
###############################################################################

temporal_degree <- function(tda_obj,...) {
	
	require(dplyr)
	require(data.table)
	require(TDAmapper)
	require(Matrix)
	
	require(foreach, quietly=T)
	require(doMC, quietly=T)
	registerDoMC(cores=20)
	
	ptl <- sum(unlist(lapply(tda_obj$points_in_level, length)))
	
	nodes <- mapperVertices(tda_obj, 1:ptl)
	adj_mat <- tda_obj$adjacency
	
	dat_list <- foreach(i=1:ncol(adj_mat)) %dopar% {
		
		sel <- which(adj_mat[,i]==1)
		
		focal_node <- apply(cbind(nodes[i,]),1, function(x) {
					str.tmp <- strsplit(as.character(x['Nodename']), split=c('V',':'))[[1]][2]
					substr(str.tmp, start=regexpr(':', str.tmp)[1], stop=regexpr(':', str.tmp)[1]) <- ','
					as.vector(as.numeric(unlist(strsplit(str.tmp,","))))[-1]
				})
		
		dat_focal <- data.frame(focal=as.vector(focal_node)) %>% group_by(focal) %>% summarise(freq=n())
		
		if(length(sel)>0) {
			
			linked_nodes <- apply(cbind(nodes[sel,]),1, function(x) {
						str.tmp <- strsplit(as.character(x['Nodename']), split=c('V',':'))[[1]][2]
						substr(str.tmp, start=regexpr(':', str.tmp)[1], stop=regexpr(':', str.tmp)[1]) <- ','
						focal_times=as.vector(as.numeric(unlist(strsplit(str.tmp,","))))[-1]
					})
			
			if(is.list(linked_nodes)) {
				dat_linked <- data.frame(focal=unlist(linked_nodes)) %>% group_by(focal) %>% summarise(freq=n())
			}
			
			else {
				dat_linked <- data.frame(focal=as.vector(linked_nodes)) %>% group_by(focal) %>% summarise(freq=n())
			} 
			
			dat_merged <- bind_rows(dat_focal,dat_linked) %>% group_by(focal) %>% summarise(freq=n()) %>% 
					filter(focal%in%c(dat_focal$focal, dat_linked$focal))
			
			mat_merged <- matrix(dat_merged$freq, nrow=nrow(dat_merged),
					ncol=nrow(dat_merged))
			diag(mat_merged) <- diag(mat_merged)-1
			
			if(length(mat_merged)==1) {
				mat_merged <- as.matrix(mat_merged)
				rownames(mat_merged) <- colnames(mat_merged) <- dat_merged$focal
			}
			
			else {
				rownames(mat_merged) <- colnames(mat_merged) <- dat_merged$focal
				mat_merged <- mat_merged[,which(as.numeric(colnames(mat_merged))%in%dat_focal$focal)]
			}
			
		}
		
		else {
			
			mat_merged <- matrix(dat_focal$freq, nrow=nrow(dat_focal),
					ncol=nrow(dat_focal))
			rownames(mat_merged) <- colnames(mat_merged) <- dat_focal$focal
			diag(mat_merged) <- diag(mat_merged)-1
		}
		
		out <- as.matrix(mat_merged)
		colnames(out) <- dat_focal$focal
		out_dat <- data.frame(focal=as.numeric(names(colSums(out))),
				freq=as.vector(colSums(out)))
		out_dat
	}
	
	
	dat_agg <- bind_rows(dat_list)
	degree_dat <- dat_agg %>% group_by(focal) %>% summarise(freq=sum(freq))
	dd <- degree_dat$freq
	# normalize
	#(dd-min(dd))/(max(dd)-min(dd))
	nms <- unique(unlist(tda_obj$points_in_vertex))
	nms <- nms[order(nms)]
	names(dd) <- nms
	dd
}

