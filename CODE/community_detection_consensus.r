#########################################################################################
# Consensus heatmap from 5 different community detection methods
#########################################################################################
# Dependencies for variables
# MAIN.r, create_discovery_set.r
#########################################################################################
mat = weightMat
ux = unlist(strsplit(rownames(mat),"\\."))
prot1 = as.numeric(ux[c(T,F)])
prot2 = as.numeric(ux[c(F,T)])
uniqs = sort(prots[unique(c(prot1,prot2))]) # unique protein names in the network

# Create graph object using the igraph library
library(igraph)
xedge = cbind(prots[prot1],prots[prot2])
g = graph_from_data_frame(xedge,directed=F)

# Community detection methods
# In decreasing order of the resulting modularity score 
# 1) cluster_spinglass (modularity:0.44, stochastic)
# 2) cluster_infomap (modularity:0434, stochastic)
# 3) cluster_louvain (modularity:0429, deterministic)
# 4) cluster_fast_greedy (modularity:0411, deterministic)
# 5) cluster_walktrap (modularity:04105, deterministic)

set.seed(111)
obj1 = cluster_spinglass(g)
obj = obj7
modularity(obj) # 0.44 
sizes(obj) # 8 modules
comm = communities(obj)
A1 = convertCommunityToAdjacency(comm,uniqs)

set.seed(111)
obj2 = cluster_infomap(g,nb.trials = 100) # stochastic
obj = obj2
modularity(obj) # 0.434 but changes every time
sizes(obj) # 11 modules
comm = communities(obj)
A2 = convertCommunityToAdjacency(comm,uniqs)

obj3 = cluster_louvain(g)
obj = obj3
modularity(obj) # 0.429 
sizes(obj) # 6 modules
comm = communities(obj)
A3 = convertCommunityToAdjacency(comm,uniqs)

obj4 = cluster_fast_greedy(g) 
obj = obj4
modularity(obj) # 0.411
sizes(obj) # 6 modules
comm = communities(obj)
A4 = convertCommunityToAdjacency(comm,uniqs)

obj5 = cluster_walktrap(g)
obj = obj5
modularity(obj) # 0.4105 
sizes(obj) # 22 modules
comm = communities(obj)
A5 = convertCommunityToAdjacency(comm,uniqs)
#####################################
# Get frequency for each pair of antibodies occurring in the same community
A = (A1 + A2 + A3 + A4 + A5) / 5
###################################################################
# Heatmap colors and linkage
white=rgb(1,1,1);red = rgb(1,0,0)
whitered = colorRampPalette(c(white,red)) 
# Auxiliary function for Ward linkage hierchical clustering
hclustWard = function(d){
  hclust(d,method="ward.D2") # this will square distances
}# end function
###################################################################
# Plot heatmap
colvec = colorblindFriendly1()[1:6]
hcx = hclust(dist(A),method="ward.D2")
classes = cutree(hcx,6)
rowAnnot = matrix(NA,nrow=1,ncol=ncol(A))
colAnnot = matrix(NA,ncol=1,nrow=nrow(A))

for(i in 1:6){
  w = which(classes==i)
  rowAnnot[w] = colvec[i]
  colAnnot[w] = colvec[i]
}#end for i

pdf("community_detection_consensus.pdf",width=9,height=9)
heatmap.3(A,hclustfun=hclustWard,col=whitered,lhei=c(0.6,4),lwid=c(0.6,4),
          cexCol=0.3,cexRow=0.3,margins=c(4,4),KeyValueName="Frequency",
          ColSideColors = colAnnot,RowSideColors = rowAnnot,side.height.fraction=0.15)
dev.off()
###################################################################









