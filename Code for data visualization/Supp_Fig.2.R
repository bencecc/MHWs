#### ---- Reproduce Supplementary Figure 2. Plots of networks ------ ####
#### ---- annotated with selected attributes. ---------------------- ####

# ---- Require libraries
require(dplyr)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)
require(TDAmapper)
require(igraph)

# ---- Load functions
source('~/MHWs/code/plot_net.R') #function to generate network plots from TDA objects
source('~/MHWs/code/color_nodes_time.R') # function to color nodes in the network

# ---- Observed SST
load('~/MHWs/data/tda_glob_res.RData') # tda data
load('~/MHWs/data/tda_glob_covar.RData') # covariate data to color nodes

# ---- Historical
load('~/MHWs/data/tda_hist_res.RData')
load('~/MHWs/data/tda_hist_covar.RData')
# ---- RCP 8.5
load('~/MHWs/data/tda_rcp85_res.RData')
load('~/MHWs/data/tda_rcp85_covar.RData')
# ---- RCP 2.6
load('~/MHWs/data/tda_rcp26_res.RData')
load('~/MHWs/data/tda_rcp26_covar.RData')

#---- Cumulative intensity
glob_intcum <- plot_net(tda_res = tda_glob_res, tda_nrows = 13514,
		col_node_dat = tda_glob_covar, col_var='intensity_cumulative', plot=F)

hist_intcum <- plot_net(tda_res = tda_hist_res, tda_nrows = 52960,
		col_node_dat = tda_hist_covar, col_var='intensity_cumulative', plot=F)

rcp85_intcum <- plot_net(tda_res = tda_rcp85_res, tda_nrows = 34675,
		col_node_dat = tda_rcp85_covar, col_var='intensity_cumulative', plot=F)

rcp26_intcum <- plot_net(tda_res = tda_rcp26_res, tda_nrows = 34675,
		col_node_dat = tda_rcp26_covar, col_var='intensity_cumulative', plot=F, is.rcp26=T)

windows(height=7, width=7)
ggarrange(glob_intcum, hist_intcum, rcp26_intcum, rcp85_intcum,
		ncol=2, nrow=2,
		common.legend=F,
		legend="none"
)

#ggsave("~/MHWs/data/cum_int_nets.pdf", useDingbats=FALSE, width = 7, height = 7)

# ---- Mean intensity
glob_meanint <- plot_net(tda_res = tda_glob_res, tda_nrows = 13514,
		col_node_dat = tda_glob_covar, col_var='intensity_mean', plot=F)

hist_meanint <- plot_net(tda_res = tda_hist_res, tda_nrows = 52960,
		col_node_dat = tda_hist_covar, col_var='intensity_mean', plot=F)

rcp85_meanint <- plot_net(tda_res = tda_rcp85_res, tda_nrows = 34675,
		col_node_dat = tda_rcp85_covar, col_var='intensity_mean', plot=F)

rcp26_meanint <- plot_net(tda_res = tda_rcp26_res, tda_nrows = 34675,
		col_node_dat = tda_rcp26_covar, col_var='intensity_mean', plot=F, is.rcp26=T)


windows(height=7, width=7)
ggarrange(glob_meanint, hist_meanint, rcp26_meanint, rcp85_meanint,
		ncol=2, nrow=2,
		#common.legend=T,
		legend="none"
)

#ggsave("~/MHWs/data/mean_int_nets.pdf", useDingbats=FALSE, width = 7, height = 7)


# ---- Number of events
glob_numevents <- plot_net(tda_res = tda_glob_res, tda_nrows = 13514,
		col_node_dat = tda_glob_covar, col_var='num_events', plot=F)

hist_numevents <- plot_net(tda_res = tda_hist_res, tda_nrows = 52960,
		col_node_dat = tda_hist_covar, col_var='num_events', plot=F)

rcp85_numevents <- plot_net(tda_res = tda_rcp85_res, tda_nrows = 34675,
		col_node_dat = tda_rcp85_covar, col_var='num_events', plot=F)

rcp26_numevents <- plot_net(tda_res = tda_rcp26_res, tda_nrows = 34675,
		col_node_dat = tda_rcp26_covar, col_var='num_events', plot=F, is.rcp26=T)


windows(height=7, width=7)
ggarrange(glob_numevents, hist_numevents, rcp26_numevents, rcp85_numevents,
		ncol=2, nrow=2,
		common.legend=F,
		legend="none"
)

#ggsave("~/MHWs/data/num_events_nets.pdf", useDingbats=FALSE, width = 7, height = 7)


		
windows(height=7, width=4)
ggarrange(glob_intcum, hist_intcum, rcp26_intcum, rcp85_intcum,
		glob_numevents, hist_numevents, rcp26_numevents, rcp85_numevents,
		ncol=2, nrow=4,
		common.legend=F,
		legend="none"
)		
		
#ggsave("~/MHWs/data/int_cum&num_events_nets.pdf", useDingbats=FALSE, width = 4, height = 7)

		
		
		
		