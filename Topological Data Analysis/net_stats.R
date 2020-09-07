# Summary stats from network. The function is called by TDA_pipeline.R
###############################################################################
net_stats <- function(net, ...) {
	
	if(class(net)!='igraph') {
		net_graph <- graph.adjacency(net$adjacency, mode="undirected")
	}
	
	else {
		net_graph <- net
	}
	
	fast_g <- fastgreedy.community(net_graph)
	fast_g_sizesComm <- sizes(fast_g)
	numComm <- length(fast_g_sizesComm)
	modularity <- modularity(fast_g)
	
	mean_degree <- mean(degree(net_graph))
	e_connectivity <- edge_connectivity(net_graph)
	
	#coreness using python
	if(class(net)=='igraph') {
		adj_graph <- as.matrix(as_adjacency_matrix(net))
	}
	
	else {
		adj_graph <- net$adjacency
	}
	
	G <- nx$to_networkx_graph(adj_graph)
	be = cp$KM_config() #cp$BE()
	be$detect(G)
	coreness = be$get_coreness()
	
	cbind(coreness=mean(unlist(coreness)), modularity=modularity, numComm=numComm, mean_deg=mean_degree, e_conn=e_connectivity)
	
}
