# December 3, 2015
# A compilation of code snippets used to obtain the results reported in Senbabaoglu et al.
# "A multi-method approach for proteomic network inference in 11 human cancers"
#########################################################################################
#setwd("~/Desktop/SanderLab/ACTIVE/NetworkDiscovery/ANALYSIS/Dec3_2015_compile_code/")
source("functions.r")
pancan = read.csv("DATA/INPUT RPPA DATA/TCGA-PANCAN11-RBN-v2_0.txt",sep="\t",header=TRUE)
types = as.character(unique(pancan[,2]))
prots = colnames(pancan)[-c(1,2)]
prots[1:8] = sub("X","",prots[1:8]) # 187 antibodies
tmp = table(pancan[,2])
sampleNums = tmp[match(types,names(tmp))] # number of samples in each tumor type

# Gold standard, returns pcInts
# pcInts is now cleaned. Self-interactions and duplications are removed. 1212 edges
load("DATA/PERA OUTPUT/bp_pancan_network_collapsed.rdata")
#########################################################################################
# Run all 13 methods with the appropriate arguments
# This code is in script run_methods.r (newSetup2.r from May5_2014_setup)
#########################################################################################
# Computing AUPR is done by the function longPrecisionRecall. For sparse methods, this
# function will randomly assign orders for the zero-weight edges
#
# The optimal parameter values are determined in the script choose_best_parameter.r
# and reported in S5 Fig.
#
# The edges and weights from the limited-recall optimal-parameter runs can be genereated
# using the print_limited_recall_edges.r script. The numbers of edges are reported in
# S6 Fig. Normalized weight for the CLR method can be obtained using the 
# normalize_CLR_weights.r script
#########################################################################################
# Principal component analysis (PCA) on the 13 network inference methods
# This step uses the M matrix (60202 edges, 13 methods), which can be generated using
# the script PCA_methods.r. The M matrix is provided as S9 Table.
#########################################################################################
# Consensus edge ranks and weights from the TOP6 methods are computed using the script
# consensus_rank_and_weight.r. These values are provided in S10 Table tabs a-k.
#########################################################################################
# The edges that have consensus edge rank less than 425 (i.e. the most significant edges)
# in each tumor type are taken and combined to form the T matrix (S11 Table)
# T is called the DISCOVERY SET and is 1008 by 11 (edges by tumor types).
# The code to obtain T is provided in the script create_discovery_set.r.
#########################################################################################
# PCA and hierarchical clustering on the tumor types are performed using the T matrix.
# This analysis is performed and the discovery set heatmap is generated with
# the script generate_discovery_set_heatmap.r
#########################################################################################
# The interactions in the discovery set are provided as input to five different 
# community detection methods. Community memberships from the output of these five
# methods are combined to arrive at a consensus matrix (S11 Table and Fig 6).
# Ward-linkage hierarchical clustering on this consensus matrix yields four robust, 
# and two less robust modules (Fig 6, S1A Table). The code can be found in community_detection_consensus.r.
#########################################################################################
# Average interaction strength for Reactome genelists, and number of matching interactions
# for each genelist are computed in the script reactome_weights_and_tracks_for_module_group.r
#########################################################################################

