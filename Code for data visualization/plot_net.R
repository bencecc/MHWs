#### ---- Function to plot TDA networks with the following arguments: ------------------------- ####
#### ---- tda_res is an object resulting from the mapper analysis (e.g., tda_glob_res); ------- ####
#### ---- tda_nrows are the rows of the dataset used by mapper (also the rows ----------------- ####
#### ---- of the col_node_dat dataframe) col_node_dat is the dataframe with ------------------- ####
#### ---- covariates annotate network nodes nodes and col_var is the selected covariate; ------ ####
#### ---- 'plot' indicates whether the network should be plotted or returned as an object. ---- ####
#### ------------------------------------------------------------------------------------------ ####

plot_net <- function(tda_res, tda_nrows, col_node_dat, col_var, plot=T,  ...) {
	
	nodes <- mapperVertices(tda_res, 1:tda_nrows)
	edges <- mapperEdges(tda_res)
	node_color_covar <- color_nodes_time(nodes_dat=nodes,
			covar_dat=col_node_dat, covar_name=col_var)
	
	if (is.rcp26) {
		set.seed(4321)
	}
	
	else {
		set.seed(1567527626)	
		}
		
	plot_graph <- graph.adjacency(tda_res$adjacency, mode="undirected")
	V(plot_graph)$color <- node_color_covar
	V(plot_graph)$node_id <- 1:nrow(nodes)
	
	my_layout <- as.data.frame(layout_nicely(plot_graph, dim = 2))
	my_layout$node_id <- V(plot_graph)$node_id
	my_layout$color <- V(plot_graph)$color
	edge_info <- get.data.frame(plot_graph) 
	edge_info$from.x <- my_layout$V1[match(edge_info$from, my_layout$node_id)]  #  match the from locations from the node data.frame we previously connected
	edge_info$from.y <- my_layout$V2[match(edge_info$from, my_layout$node_id)]
	edge_info$to.x <- my_layout$V1[match(edge_info$to, my_layout$node_id)]  #  match the to locations from the node data.frame we previously connected
	edge_info$to.y <- my_layout$V2[match(edge_info$to, my_layout$node_id)]
	
	qn <- quantile(node_color_covar, c(0.99), na.rm = TRUE)	
	col_pal <- rev(colorRampPalette(brewer.pal(11, 'Spectral'))(14))
	
	plot_obj <- ggplot(data=my_layout, aes(x=V1, y=V2))+
			geom_segment(data=edge_info, aes(x=from.x,xend = to.x, y=from.y,yend = to.y), size=0.3, colour="darkgrey")+
			geom_point(data=my_layout %>% filter(color<=qn), aes(color=color),
					alpha=1, size=0.5)+
			geom_point(data=my_layout %>% filter(color>qn), color=col_pal[14],
					alpha=1, size=0.5)+
			scale_color_gradientn(name=col_var,
					colors=rev(colorRampPalette(brewer.pal(11, 'Spectral'))(14)),
					na.value=NA, labels=scales::scientific)+
			theme_void()+
			theme(legend.text=element_text(size=11),
					legend.title=element_text(size=12),
					legend.position = 'right',
					legend.direction = 'vertical',
					legend.key.size = unit(0.6, "cm"))
	
	if(plot) {
		plot(plot_obj)
	}
	
	else plot_obj
}		
