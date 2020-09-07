# MHWs
Code to reproduce the results in: "Complex networks of marine heatwaves reveal tipping points in the global ocean".

This document describes how to reproduce the analyses presented in ‘Complex networks of marine heatwaves reveal tipping points in the global ocean’. To reproduce the results, the provided code and data should be placed in folders ~/MHWs/code and ~/MHWs/data, respectively. Data are available at https://figshare.com/s/ec9061c449031aa2b20e

The original analysis was performed using R-3.5.1 on a CentOS 7 cluster with 7 nodes and 72 cores per node. R-4.0.1 is used for data visualization. For proper parallelization, .sh scripts calling R functions and .txt files with parameter values are used. Remember to make .sh scripts executable using chmod +x filename.sh. Ensures that the directory paths and the .sh and .txt files provided are correctly uploaded on your system.

Several functions use multithreading within nodes and this is done with the ‘foreach’ package, which requires registering the number of cores. This number varies among the different functions depending on how much memory is used in performing the calculations. You may wish to consider changing this parameter or use other parallelization approaches that may be most appropriate for your computational environment.

The original CMIP5 data and the contemporary satellite SST observations are not provided here, but they are freely available at https://esgf-node.llnl.gov/search/cmip5/
and https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html, respectively. CMIP5 data were initially processed with Climate Data Operators (CDO) to generate ncdf4 files that were later converted into raster stacks in R. The scripts provided here assume these ncdf4 files have been properly collated and are available to the user to perform the full analysis of MHWs. To facilitate reproducing the results of the paper, the output files of the ‘heatwaveR’ package including all the MHWs statistics are provided. These files provide the input to perform Topological Data Analysis (TDA) and Event Synchronization (ES) and to obtain the associated statistical and visualization results. The sst_0.nc file of contemporary SSTs (1981-2018) is provided to enable running the full pipeline on remotely-sensed MHWs.

The key scripts to reproduce the pipelines are TDA_pipeline and ES_pipeline. These scripts use several R packages and call many functions in R. In addition, some of the analyses use .cpp or .py functions. The ‘Rcpp’ and ‘reticulate’ packages are required allow the use of C++ and Python resources in R. In addition, the Python modules ‘cpalgorithm’ and ‘networkx’ are needed.

Code for Topological Data Analysis (TDA)

TDA_pipeline.R – main function to perform all steps in the TDA analysis
mapper2D_parallel.R – function to perform the key steps of the TDA analysis in parallel 
eucDist.cpp – fast computation of Euclidean distances (called in mapper2D_parallel.R)
HeatWaves.R – detection of MHWs using a climatology from the same input .nc file
HeatWavesClim.R – detection of MHWs using custom data input for climatology (e.g. historical)
covar_gen.R – generates covariates to annotate networks
tda_optim.R – TDA analysis on a grid of resolution and gain parameters
net_stats.R – derive network statistics from the adjacency matrix of a TDA object (output of tda_optim)
temporal_degree.R – compute the node degree from tda object
rand_phase_net.R – null model analysis (constant or random phase randomization)
cluster_cut_off.R – function involved in hierarchical clustering
temporal_connectivity.R – function to compute the temporal connectivity matrix and node degree from a TDA object

Additional code
Each of the .R functions above has a corresponding .sh script to enable parallel computation.

Code for Event Synchronization Analysis (ES)

EvSync_pipeline.R –main function to perform all steps in the ES analysis
EvSync_tlen.R – compute the total number of events between pairs of cells/pixels
EvSync_null_model.R – set thresholds to identify significant connections through null models
EvSync_nopar.cpp – compute ES from two vectors of indices of ES occurrences (called in EvSync_null_model.R)
EvSync_main-R – compute ES between all pairs of cells and retain significant connections
EvSync_sign_serial.cpp and EvSync_clust.cpp – compute ES in parallel (called in EvSync_main.R)
EvSync_cell_dist.R – compute great-circle distances between significantly linked cells
EvSync_allDist.R – compute all possible great-circle distances

Additional code
Each of the .R functions above has a corresponding .sh script to enable parallel computation.

Code for visualization

Fig.1.R – Fig. 1 of main paper
Fig. 2.R – Fig. 2 of main paper
Fig.3.R – Fig. 3 of main paper
Fig.4.R – Fig. 4 of main paper
ExtData_Fig.2.R – Extended Data Fig. 2
ExtData_Fig.3.R – Extended Data Fig. 3
ExtData_Fig.4.R – Extended Data Fig. 4
ExtData_Fig.5.R – Extended Data Fig. 5
ExtData_Figs_6_9.R – Extended Data Figs. 6-9
ExtData_Fig.10.R – Extended Data Fig. 10
color_nodes_time.R – aggregates (mean, sum) covariate values among timeframes in a node; used to annotate network graphs by covariate values
net_to_map.R – generates .nc files of covariates to produce maps from selected nodes in a network
net_to_time.R – displays time series of node attributes from selected nodes in a network
plot_time – plot a TDA network

NOTE: Extended Data Fig. 1 is artwork produced with Adobe Illustrator

Data files (available at https://figshare.com/s/ec9061c449031aa2b20e)

The analysis requires several data files. These are divided in four groups identified by a unique prefix: glob, hist, recp26, rcp85, to reflect datasets originating from remotely-sensed, historical RCP 2.6 and RCP 8.5 SSTs, respectively. For each group, different datasets are provided, which allow to reproduce the TDA and Event Synchronization analysis and to generate the main figures in the paper. A brief description of each dataset is provided here using glob data as an example (suffix .RData is omitted from the description below for clarity).

hw_glob_df – data frame originating as the output of functions in the heatwaveR package 

tda_glob - a data table of the presence/absence of MHW occurrences (rows are days and columns are pixels)

tda_glob_covar – data frame of covariates to annotate MHW networks

tda_glob_res – output of the TDA Mapper algorithm including the adjacency matrix to build the MHW network
glob_res_mat – list including the temporal connectivity matrix and the normalized degree

glob_net_perturb_dat – list of data to visualize the 16 networks used in Extended Data Fig. 6a

tda_glob_xx_yy - output of the TDA Mapper algorithm for the extreme resolution (xx) and gain (yy) parameters used for the sensitivity analysis presented in Extended Data Fig. 6b,c

tda_glob_null_stats – list of results from random phase null models

tda_glob_nulldep_stats – list of results from constant phase null models

glob_optimxxyy_null_stats – list of results from null models on extreme networks identified by their resolution (xx) and gain (yy) parameters
 
glob_years_lens – tibble including elapsed days for plotting

tcm_glob_red – temporal connectivity matrix (data table) reduced by a factor 10 for plotting

Additional datasets for analysis

Datasets of node degree and results from null distributions are provided for the Historical and RCP 8.5 scenarios to reproduce Extended Data Figs. 3 and 4:

hist_degree_ind/ rcp85_degree_ind – node degree from individual models with mean and SD from null distributions to generate CIs.

hist_deg_null/rcp85_deg_null – mean and SD from null models to compute CIs of node degree

Additional data for visualization

glob_map_duration.nc – this .nc file includes raster layers of average duration of MHWs for each node in the corresponding network (analogous files are available for hist, rcp26 and rcp85) and allow quick coloring of maps. Please be aware that the Time axis identifies network nodes, not time in the .nc files. For example, the rcp85 network has 1922 nodes, so the corresponding rcp85_map_duration.nc file has 1922 layers. The file can be red into a raster brick (or stack) to reproduce Extended Data Fig. 2. Only one of eight possible covariates is provided to save disk space. Using function net_to_map it is possible to generate .nc files for the whole set of covariates. This is also described in the ExtData_Fig.2.R script. These .nc files allow to quickly select mapping variables and are used in the web application available at: http://calcoloecologia.biologia.unipi.it:3838/MHW_App.

es_plotdat_tauXX – datasets to plot the results of event synchronization analysis

spread_year_hist –dataset of the spread of years along the temporal connectivity matrix columns used to generate Fig. 3.

world_dt_robin – data table with coastlines for global plots of MHWs with robin projection 

Required R libraries

abind
bigstatsr
cmocean
colorspace
data.table
doMC
dplyr
fastcluster
foreach
fractal
ggpubr
ggplot2
ggsci
ggthemes
heatwaveR
igraph
maps
Matrix
ncdf4
ncdf4.helpers
PCICt
raster
rasterVis
RColorBrewer
Rcpp
rgdal
rgeos
reticulate
sf
stars
stringr
TDAmapper
tidyr
uwot
VertexSort
viridis
wesanderson
wordspace

