# Derivation of temporal connectivity matrix; the function uses a tda_res
# file obtained as an output of a TDA analysis (e.g., tda_optim).
###############################################################################

temporal_connectivity <- function(tda_obj,...) {
	
	require(dplyr)
	require(data.table)
	require(TDAmapper)
	require(Matrix)
	require(foreach, quietly=T)
	require(doMC, quietly=T)
	registerDoMC(cores=30)
	
	ptl <- sum(unlist(lapply(tda_obj$points_in_level, length)))
	
	nodes <- mapperVertices(tda_obj, 1:ptl)
	adj_mat <- tda_obj$adjacency
	
	mat_list <- foreach(i=1:ncol(adj_mat)) %dopar% {
		
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
		out
	}
	
	library(doParallel)
	cl <- makeCluster(30)
	chunks <- clusterSplit(cl, 1:length(mat_list))
	ssp <- lapply(chunks, function(x) mat_list[x])
	stopCluster(cl)
	
	tmp_mat <- foreach(i=1:length(ssp)) %dopar% {
		res <- rbindlist(lapply(ssp[[i]], function(x) data.table::setDT(as.data.frame(x), 
									keep.rownames = TRUE)), fill = TRUE, use.names=T)[, lapply(.SD, sum, na.rm = TRUE), by = rn]
	}
	
	tmp_mat <- rbindlist(tmp_mat, fill=T)     
	
	setkey(tmp_mat, "rn")
	tmp_zero <- setDT(tmp_mat)[,lapply(.SD, function(x) ifelse(is.na(x),0,x))][,lapply(.SD, sum), by='rn']

	row_order <- sort(as.numeric(tmp_zero$rn))
	row_idx <- match(row_order, as.numeric(tmp_zero$rn))
	tmp_zero <- tmp_zero[row_idx,]

	tmp_zero <- tmp_zero[,-'rn']
	
	col_order <- sort(as.numeric(colnames(tmp_zero)))
	col_idx <- match(col_order,as.numeric(colnames(tmp_zero)))
	tmp_ord <- data.table::setcolorder(tmp_zero, col_idx)
	deg_tmp <- apply(tmp_ord,2,function(x) sum(x))
	degree_dat <- (deg_tmp-min(deg_tmp))/(max(deg_tmp)-min(deg_tmp))
	
	tmp_ord <- Matrix(as.matrix(tmp_ord))
	colnames(tmp_ord) <- colnames(tmp_zero)
	tcm_mat <- (tmp_ord-min(tmp_ord))/(max(tmp_ord)-min(tmp_ord))
	res_mat <- list(tcm_mat=tcm_mat, norm_degree=degree_dat)
	res_mat
}

#system.time(glob_res_mat <- temporal_connectivity(tda_glob_res))
#save(glob_res_mat, file='glob_res_mat.RData')



