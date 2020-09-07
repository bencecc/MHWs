# This function decides where to cut the hierarchical clustering tree to
# define clusters within a level set. This is part of Paul Pearson's et al
# TDA Mapper package (https://github.com/paultpearson/TDAmapper).
# The function is called by mapper2D_parallel
###############################################################################

cluster_cutoff_at_first_empty_bin <- function(heights, diam, num_bins_when_clustering) {
	
	# if there are only two points (one height value), then we have a single cluster
	if (length(heights) == 1) {
		if (heights == diam) {
			cutoff <- Inf
			return(cutoff)
		}
	}
	
	bin_breaks <- seq(from=min(heights), to=diam, 
			by=(diam - min(heights))/num_bins_when_clustering)
	if (length(bin_breaks) == 1) { bin_breaks <- 1 }
	
	myhist <- hist(c(heights,diam), breaks=bin_breaks, plot=FALSE)
	z <- (myhist$counts == 0)
	if (sum(z) == 0) {
		cutoff <- Inf
		return(cutoff)
	} else {
		#  which returns the indices of the logical vector (z == TRUE), min gives the smallest index
		cutoff <- myhist$mids[ min(which(z == TRUE)) ]
		return(cutoff)
	}
	
}

