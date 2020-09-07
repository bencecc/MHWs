# This is a parallelized version of Paul Pearson's et al TDA Mapper function;
# (https://github.com/paultpearson/TDAmapper). It uses a .cpp function for fast
# calculation of Euclidean distances. The function is called by rand_phase_net.R,
# and tda_optim.R in TDA_pipeline.R. The output includes the adjacency matrix
# used to build networks and to compute network statistics.
#################################################################################
mapper2D_parallel <- function(
		data,
		is.distance=F,
		ncores=NULL,
		filter_values = NULL,
		num_intervals = c(5,5),
		percent_overlap = 50,
		num_bins_when_clustering = 10
) {
	
	require(foreach, quietly=T)
	require(doMC, quietly=T)
	registerDoMC(cores=ncores)
	require(wordspace)
	require(data.table)
	require(Rclusterpp)
	require(sp)
	require(Rcpp)
	
	source('~/MHWs/code/cluster_cut_off.R')
	sourceCpp('~/MHWs/code/eucDist.cpp')
	
	######################################################################
	# initialize variables
	vertex_index <- 0
	
	# indexed from 1 to the number of vertices
	level_of_vertex <- c()
	points_in_vertex <- list()
	
	vertices_in_level <- list()
	
	filter_min_1 <- min(filter_values[[1]])
	filter_max_1 <- max(filter_values[[1]])
	filter_min_2 <- min(filter_values[[2]])
	filter_max_2 <- max(filter_values[[2]])
	
	interval_length_1 <- (filter_max_1 - filter_min_1) / (num_intervals[1] - (num_intervals[1] - 1) * percent_overlap/100 )
	interval_length_2 <- (filter_max_2 - filter_min_2) / (num_intervals[2] - (num_intervals[2] - 1) * percent_overlap/100 )
	
	step_size_1 <- interval_length_1 * (1 - percent_overlap/100)
	step_size_2 <- interval_length_2 * (1 - percent_overlap/100)
	
	num_levels <- num_intervals[1] * num_intervals[2]
	
	level_indices_1 <- rep(1:num_intervals[1], num_intervals[2])
	level_indices_2 <- rep(1:num_intervals[2], each=num_intervals[1])
	
	# begin mapper main loop
	out_ini <- foreach(level=1:num_levels) %do%  {
		
		level_1 <- level_indices_1[level]
		level_2 <- level_indices_2[level]
		
		min_value_in_level_1 <- filter_min_1 + (level_1 - 1) * step_size_1
		min_value_in_level_2 <- filter_min_2 + (level_2 - 1) * step_size_2
		max_value_in_level_1 <- min_value_in_level_1 + interval_length_1
		max_value_in_level_2 <- min_value_in_level_2 + interval_length_2
		
		points_in_level_logical <-
				(min_value_in_level_1 <= filter_values[[1]]) &
				(min_value_in_level_2 <= filter_values[[2]]) &
				(filter_values[[1]] <= max_value_in_level_1) &
				(filter_values[[2]] <= max_value_in_level_2)
		
		num_points_in_level <- sum(points_in_level_logical)
		points_in_level <- which(points_in_level_logical==TRUE)
		
		list(points_in_level_logical=points_in_level_logical, num_points_in_level=num_points_in_level,
				points_in_level=points_in_level)
	}
	
	out <- foreach(i=1:length(out_ini)) %dopar% {
		
		points_in_level_logical <- out_ini[[i]][[1]]
		num_points_in_level <- out_ini[[i]][[2]]
		points_in_level <- out_ini[[i]][[3]]
		
		if (num_points_in_level == 0) {
			num_vertices_in_level <- -1
			cluster_indices_within_level <- -1
		}
		
		else {
			
			if (num_points_in_level == 1) {
				num_vertices_in_level <- 1
				cluster_indices_within_level <- 1
				
			}
			
			if (num_points_in_level > 1) {
				
				if(!is.distance){
					
					# using eucDist_parallel to obtain euclidean distances
					mat <- as.matrix(data[points_in_level_logical])
					level_distance_matrix <- as.dist(eucDist_parallel(mat))
					rm(mat)
					
				}
				
				else {
					level_distance_matrix <- Matrix(as.matrix(data)[points_in_level_logical,points_in_level_logical])
				}
				
				level_max_distance <- max(level_distance_matrix)
				
				level_hcluster_ouput <- Rclusterpp.hclust(level_distance_matrix, method="single")
				
				heights <- level_hcluster_ouput$height
				cutoff <- cluster_cutoff_at_first_empty_bin(heights, level_max_distance, num_bins_when_clustering)
				
				cluster_indices_within_level <- as.vector(cutree(level_hcluster_ouput, h=cutoff))
				num_vertices_in_level <- max(cluster_indices_within_level)
			}
		}
		
		res <- list(num_vertices_in_level=num_vertices_in_level,
				points_in_level=points_in_level, points_in_level_logical=points_in_level_logical,
				cluster_indices_within_level=cluster_indices_within_level)
		res
		
	}	
	
	
	level <- 1:length(out)
	
	num_vertices_in_level <- lapply(out, function(x) x[[1]])
	points_in_level <- sapply(out, function(x) x[[2]])
	points_in_level_logical <- lapply(out, function(x) x[[3]])
	cluster_indices_within_level <- lapply(out, function(x) x[[4]])
	
	for (i in 1:length(level)) {
		
		if(num_vertices_in_level[[i]]==-1||is.null(num_vertices_in_level[[i]])) {
			vertices_in_level[[i]] <- -1
			points_in_level[[i]] <- as.integer(0)
			cluster_indices_within_level[[i]] <- -1
		}	
		
		else {
			vertices_in_level[[i]] <- vertex_index + 1:num_vertices_in_level[[i]]
			
			for (j in 1:num_vertices_in_level[[i]]) {
				
				vertex_index <- vertex_index + 1
				
				level_of_vertex[vertex_index] <- level[i]
				
				points_in_vertex[[vertex_index]] <- which(points_in_level_logical[[i]]==TRUE)[cluster_indices_within_level[[i]] == j]
				
			}
		}
		
	}
# end mapper main loop
	
	# Note: num_vertices = vertex index.
	# Create the adjacency matrix for the graph, starting with a matrix of zeros
	adja <- mat.or.vec( vertex_index, vertex_index )
	
	for (i in 1:num_intervals[1]) {
		for (j in 2:num_intervals[2]) {
			
			k1 <- which( (level_indices_1 == i) & (level_indices_2 == j) )
			k2 <- which( (level_indices_1 == i) & (level_indices_2 == j-1))
			
			if ( (vertices_in_level[[k1]][1] != -1) & (vertices_in_level[[k2]][1] != -1) ) {
				
				for (v1 in vertices_in_level[[k1]]) {
					for (v2 in vertices_in_level[[k2]]) {
						adja[v1,v2] <- ( length(intersect(points_in_vertex[[v1]],
													points_in_vertex[[v2]])) > 0 )
						adja[v2,v1] <- adja[v1,v2]
					}
				}
				
			}
		} 
	}
	
	for (j in 1:num_intervals[2]) {
		for (i in 2:num_intervals[1]) {
			
			k1 <- which( (level_indices_1 == i) & (level_indices_2 == j) )
			k2 <- which( (level_indices_1 == i-1) & (level_indices_2 == j))
			
			if ( (vertices_in_level[[k1]][1] != -1) & (vertices_in_level[[k2]][1] != -1) ) 	{
				for (v1 in vertices_in_level[[k1]]) {
					for (v2 in vertices_in_level[[k2]]) {
						adja[v1,v2] <- ( length(intersect(points_in_vertex[[v1]],
													points_in_vertex[[v2]])) > 0 )
						adja[v2,v1] <- adja[v1,v2]
					}
				}
				
			}
		} 
	}
	
	mapperoutput <- list(adjacency = adja,
			num_vertices = vertex_index,
			level_of_vertex = level_of_vertex,
			points_in_vertex = points_in_vertex,
			points_in_level = points_in_level,
			vertices_in_level = vertices_in_level
	)
	
	class(mapperoutput) <- "TDAmapper"
	
	return(mapperoutput)
	
}



