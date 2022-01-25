##################################################################################################
# Choosing the best parameter value for methods that require user optimization
##################################################################################################
rm(list=ls())
homedir <- "/Users/senbabay/Documents/AKORKUT"
setwd(file.path(homedir, "2022_01_16_brca_intro/"))

source(file.path(homedir,"CODE/functions.r"))
load(file.path(homedir,"DATA/PERA/bp_pancan_network_collapsed.rdata")) # ground truth from PERA
brca <- read.csv(file.path(homedir,"DATA/TCGA-BRCA-L4_protnet.txt"),sep="\t",header=TRUE)
#################################################
# The column name for sample identifiers
sampleID <- "SampleID"
# Antibody names
prots <- setdiff(colnames(brca),sampleID)
# If antibody names start with X to avoid a number as initial character, remove X
starts.with.x <- grep("^X[0-9]",prots)
prots[starts.with.x] <- sub("X","",prots[starts.with.x])

##################################################################################################
# GLASSO
##################################################################################################
edgelist <- readRDS(file.path(homedir,"OUTPUT/EDGELISTS/glasso.rds"))
areaVec <- rep(NA,length(edgelist))
sortedEdges <- vector("list",length(edgelist))

for(j in 1:length(edgelist)){
  if(j %% 20 == 0){
    print(j)
  }
  res <- longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE)
  areaVec[j] <- res$area
  sortedEdges[[j]] <- res$sortedEdges
}#end for j
bestRhoInd <- which.max(areaVec)
maxArea <- max(areaVec)
bestEdges <- sortedEdges[[bestRhoInd]]

# edges with significance
edges <- edgelist[[bestRhoInd]]
edges$node1name <- prots[edges$node1]
edges$node2name <- prots[edges$node2]

write.table(edges,file.path(homedir,"OUTPUT/predictions_glasso.csv"),sep=",",quote=FALSE,row.names=FALSE)

##################################################################################################
# PLSNET
##################################################################################################
edgelist <- readRDS(file.path(homedir,"OUTPUT/EDGELISTS/plsnet.rds"))

areaVec <- rep(NA,length(edgelist))
sortedEdges <- vector("list",length(edgelist))
kvec <- as.numeric(gsub("k","",names(edgelist)))

for(j in 1:length(edgelist)){
  res <- longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE)
  areaVec[j] <- res$area
  sortedEdges[[j]] <- res$sortedEdges
}#end for j

bestKInd <- which.max(areaVec)
bestK <- kvec[bestKInd]
maxArea <- max(areaVec)
bestEdges <- sortedEdges[[bestKInd]]

# edges with significance
edges <- edgelist[[bestKInd]]
edges$node1name <- prots[edges$node1]
edges$node2name <- prots[edges$node2]

write.table(edges,file.path(homedir,"OUTPUT/predictions_plsnet.csv"),sep=",",quote=FALSE,row.names=FALSE)
##################################################################################################
# RIDGENET
##################################################################################################
edgelist <- readRDS(file.path(homedir,"OUTPUT/EDGELISTS/ridgenet.rds"))

areaVec <- rep(NA,length(edgelist))
sortedEdges <- vector("list",length(edgelist))
kvec <- as.numeric(gsub("k","",names(edgelist)))

for(j in 1:length(edgelist)){
  res <- longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE)
  areaVec[j] <- res$area
  sortedEdges[[j]] <- res$sortedEdges
}#end for j

bestKInd <- which.max(areaVec)
bestK <- kvec[bestKInd]
maxArea <- max(areaVec)
bestEdges <- sortedEdges[[bestKInd]]

# edges with significance
edges <- edgelist[[bestKInd]]
edges$node1name <- prots[edges$node1]
edges$node2name <- prots[edges$node2]

write.table(edges,file.path(homedir,"OUTPUT/predictions_ridgenet.csv"),sep=",",quote=FALSE,row.names=FALSE)

##################################################################################################
# LASSONET
##################################################################################################
edgelist <- readRDS(file.path(homedir,"OUTPUT/EDGELISTS/lassonet.rds"))
areaVec <- rep(NA,length(edgelist))
sortedEdges <- vector("list",length(edgelist))
lambdaVec <- as.numeric(gsub("lambda","",names(edgelist)))

for(j in 1:length(edgelist)){
  if(j %% 20 == 0){
    print(j)
  }
  res <- longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE)
  areaVec[j] <- res$area
  sortedEdges[[j]] <- res$sortedEdges
}#end for j

bestLInd <- which.max(areaVec)
bestL <- kvec[bestLInd]
maxArea <- max(areaVec)
bestEdges <- sortedEdges[[bestLInd]]

# edges with significance
edges <- edgelist[[bestLInd]]
edges$node1name <- prots[edges$node1]
edges$node2name <- prots[edges$node2]

write.table(edges,file.path(homedir,"OUTPUT/predictions_lassonet.csv"),sep=",",quote=FALSE,row.names=FALSE)

##################################################################################################
# ARACNE.A
##################################################################################################
edgelist <- readRDS(file.path(homedir,"OUTPUT/EDGELISTS/aracne_a.rds"))
knnvec <- as.numeric(gsub("k","",names(edgelist)))
epsvec <- as.numeric(gsub("eps","",names(edgelist[[1]])))

# Optimize k and epsilon together (grid search)
# meaning find the maximum area in the matrix of k and epsilon
maxArea <- 0
for(j in 1:length(knnvec)){
  print(knnvec[j])
  edgefile <- edgelist[[j]]
  for(w in 1:length(edgefile)){
    if(w %% 20 == 0){
      print(w)
    }
    res <- longPrecisionRecall(edgefile[[w]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE,node1ind=2,node2ind=3)
    if(res$area > maxArea){
      maxArea <- res$area
      bestArea <- res$area
      bestK <- knnvec[j]
      bestEps <- epsvec[w]
      bestEdges <- res$sortedEdges
      
      bestk_ix <- j
      besteps_ix <- w
    }#end if
  }#end for w
}#end for j

# edges with significance
edges <- as.data.frame(edgelist[[bestk_ix]][[besteps_ix]])
edges$node1name <- prots[edges$row]
edges$node2name <- prots[edges$col]

write.table(edges,file.path(homedir,"OUTPUT/predictions_aracne_additive.csv"),sep=",",quote=FALSE,row.names=FALSE)

##################################################################################################
# ARACNE.M
##################################################################################################
edgelist <- readRDS(file.path(homedir,"OUTPUT/EDGELISTS/aracne_m.rds"))

knnvec <- as.numeric(gsub("k","",names(edgelist)))
tauvec <- as.numeric(gsub("tau","",names(edgelist[[1]])))

# Optimize k and tau together (grid search)
# meaning find the maximum area in the matrix of k and tau
maxArea<-0
for(j in 1:length(knnvec)){
  print(knnvec[j])
  edgefile<-edgelist[[j]]
  for(w in 1:length(edgefile)){
    if(w %% 20 == 0){
      print(w)
    }
    res<-longPrecisionRecall(edgefile[[w]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE,node1ind=2,node2ind=3)
    if(res$area > maxArea){
      maxArea<-res$area
      bestArea<-res$area
      bestK<-knnvec[j]
      bestTau<-tauvec[w]
      bestEdges<-res$sortedEdges
      
      bestk_ix <- j
      besttau_ix <- w
    }#end if
  }#end for w
}#end for j

# edges with significance
edges <- as.data.frame(edgelist[[bestk_ix]][[besttau_ix]])
edges$node1name <- prots[edges$row]
edges$node2name <- prots[edges$col]

write.table(edges,file.path(homedir,"OUTPUT/predictions_aracne_multiplicative.csv"),sep=",",quote=FALSE,row.names=FALSE)

##################################################################################################
# CLR
##################################################################################################
edgelist <- readRDS(file.path(homedir,"OUTPUT/EDGELISTS/CLR.rds"))
areaVec <- rep(NA,length(edgelist))
sortedEdges <- vector("list",length(edgelist))
knnvec <- as.numeric(gsub("k","",names(edgelist)))

for(j in 1:length(edgelist)){
  res <- longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE)
  areaVec[j] <- res$area
  sortedEdges[[j]] <- res$sortedEdges
}#end for j

bestKInd <- which.max(areaVec)
bestK <- knnvec[bestKInd]
maxArea <- max(areaVec)
bestEdges <- sortedEdges[[bestKInd]]

# edges with significance
edges <- as.data.frame(edgelist[[bestKInd]])
edges$node1name <- prots[edges$row]
edges$node2name <- prots[edges$col]

write.table(edges,file.path(homedir,"OUTPUT/predictions_CLR.csv"),sep=",",quote=FALSE,row.names=FALSE)

##################################################################################################
# MRNET
##################################################################################################
edgelist <- readRDS(file.path(homedir,"OUTPUT/EDGELISTS/MRNET.rds"))
areaVec <- rep(NA,length(edgelist))
sortedEdges <- vector("list",length(edgelist))
knnvec <- as.numeric(gsub("k","",names(edgelist)))

for(j in 1:length(edgelist)){
  res <- longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE)
  areaVec[j] <- res$area
  sortedEdges[[j]] <- res$sortedEdges
}#end for j

bestKInd <- which.max(areaVec)
bestK <- knnvec[bestKInd]
maxArea <- max(areaVec)
bestEdges <- sortedEdges[[bestKInd]]

# edges with significance
edges <- as.data.frame(edgelist[[bestKInd]])
edges$node1name <- prots[edges$row]
edges$node2name <- prots[edges$col]

write.table(edges,file.path(homedir,"OUTPUT/predictions_MRNET.csv"),sep=",",quote=FALSE,row.names=FALSE)
