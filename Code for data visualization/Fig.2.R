#### ---- Generates networks to compose Fig. 2 ---- ####

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
		col_node_dat = tda_glob_covar, col_var='duration', plot=F)

hist_intcum <- plot_net(tda_res = tda_hist_res, tda_nrows = 52960,
		col_node_dat = tda_hist_covar, col_var='duration', plot=F)

rcp85_intcum <- plot_net(tda_res = tda_rcp85_res, tda_nrows = 34675,
		col_node_dat = tda_rcp85_covar, col_var='duration', plot=F)

rcp26_intcum <- plot_net(tda_res = tda_rcp26_res, tda_nrows = 34675,
		col_node_dat = tda_rcp26_covar, col_var='duration', plot=F, is.rcp26=T)

windows(height=7, width=7)
ggarrange(glob_intcum, hist_intcum, rcp26_intcum, rcp85_intcum,
		ncol=2, nrow=2,
		common.legend=F,
		legend="none"
)

#ggsave("~/MHWs/data/Fig2_nets.pdf", useDingbats=FALSE, width = 7, height = 7)

