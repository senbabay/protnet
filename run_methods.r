#################################################
# Run 13 network inference methods
#################################################
# Required libraries and functions
#################################################
# ppls and parcor are removed from CRAN, they need to be downloaded and installed
#install.packages("ppls_1.6-1.tar.gz",repos=NULL,type="source",dependencies=TRUE)
#install.packages("parcor_0.2-6.tar.gz",repos=NULL,type="source",dependencies=TRUE)
library(parcor)
library(glasso)
library(parmigene) # mutual information methods
library(magrittr)
library(dplyr)
rm(list=ls())
#################################################
# Define paths and read in data
#################################################
homedir <- "/Users/senbabay/Documents/AKORKUT"
setwd(file.path(homedir, "2022_01_16_brca_intro/"))
dir.create("OUTPUT")
dir.create("OUTPUT/EDGELISTS")

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
# Create data matrix from input data (some methods only work with matrix input)
rppa <- brca %>% dplyr::select(-matches(sampleID)) %>% as.matrix(.)
rownames(rppa) <- brca[[sampleID]]
#################################################
# Methods that do not require parameter optimization
#################################################
## * GeneNet (shrunk inverse covariance)
## * Spearman correlation
## * Pearson correlation
## * Simple partial correlation (regular inverse covariance)
#################################################
# GeneNet
#################################################
res <- ggm.estimate.pcor(rppa)
edges = network.test.edges(res,verbose=FALSE,plot=FALSE)
edges$node1name <- prots[edges$node1]
edges$node2name <- prots[edges$node2]
write.table(edges,file.path(homedir,"OUTPUT/predictions_genenet.csv"),sep=",",quote=FALSE,row.names=FALSE)

#################################################
### Spearman correlation
#################################################
res <- cor(rppa,method="spearman")
edges = network.test.edges(res,verbose=FALSE,plot=FALSE)
edges$node1name <- prots[edges$node1]
edges$node2name <- prots[edges$node2]
write.table(edges,file.path(homedir,"OUTPUT/predictions_spearmanCor.csv"),sep=",",quote=FALSE,row.names=FALSE)

#################################################
# Pearson correlation
#################################################
res <- cor(rppa,method="pearson")
edges = network.test.edges(res,verbose=FALSE,plot=FALSE)
edges$node1name <- prots[edges$node1]
edges$node2name <- prots[edges$node2]
write.table(edges,file.path(homedir,"OUTPUT/predictions_pearsonCor.csv"),sep=",",quote=FALSE,row.names=FALSE)

#################################################
# Simple partial correlation
#################################################
res <- cor2pcor(cor(rppa,method="pearson"))
edges = network.test.edges(res,verbose=FALSE,plot=FALSE)
edges$node1name <- prots[edges$node1]
edges$node2name <- prots[edges$node2]
write.table(edges,file.path(homedir,"OUTPUT/predictions_simplePartialCorrelation.csv"),sep=",",quote=FALSE,row.names=FALSE)

#################################################
# Methods that require parameter optimization
#################################################
## * GLASSO (graphical lasso)
## * PLSNET (Partial least squares)
## * RidgeNET
## * LassoNET
#################################################
### GLASSO
#################################################
network.method = "glasso"

# GLASSO auxiliary function
runGlasso = function(covrppa,rho=0.01){
  gl <- glasso(covrppa,rho=rho)
  res <- invcov2parcor(gl$wi)
  return(res)
}#end function

# Vary parameter rho
# 0 not included, causes convergence problems
rhoVec <- setdiff(seq(0,1,length.out=150),0)

### Read in RPPA files and get p-values for all edges
edgeset <- vector("list",length(rhoVec))
covrppa <- cov(rppa)

for(j in 1:length(rhoVec)){
  print(j)
  # Get matrix of weights
  res <- runGlasso(covrppa,rho=rhoVec[j])
  edges = network.test.edges(res,verbose=FALSE,plot=FALSE)
  edgeset[[j]] = edges[which(abs(edges[,1]) > 0),]
}#end for j

# Remove list elements that have 0 rows (no edges predicted)
names(edgeset) <- paste0("rho",rhoVec)
keep.ix <- which(lapply(edgeset,nrow) > 0)
edgeset <- edgeset[keep.ix]
saveRDS(edgeset,paste0(homedir,"/OUTPUT/EDGELISTS",network.method,".rds"))

#################################################
### Partial Least Squares
#################################################
# Vary parameter K, keep ncomp fixed at 30
kvec <- c(3,4,5,10,20)
ncomp <- 30
network.method <- "plsnet"

edgeset = vector("list",length(kvec))
for(j in 1:length(kvec)){
  cat("j is",j,"\n")
  plsfit <- pls.net(rppa,scale=TRUE,k=kvec[j],ncomp=ncomp,verbose=FALSE)
  edgeset[[j]] = network.test.edges(plsfit$pcor,verbose=FALSE,plot=FALSE)
  print(plsfit$m)
}#end for j

# Remove list elements that have 0 rows (no edges predicted)
names(edgeset) <- paste0("k",kvec)
keep.ix <- which(lapply(edgeset,nrow) > 0)
edgeset <- edgeset[keep.ix]
saveRDS(edgeset,paste0(homedir,"/OUTPUT/EDGELISTS",network.method,".rds"))

#################################################
### RIDGENET
#################################################
# Vary parameter K
kvec = c(3,4,5,10,20)
network.method = "ridgenet"

edgeset = vector("list",length(kvec))
for(j in 1:length(kvec)){
  cat("j is",j,"\n")
  ridgefit <- ridge.net.vLambda(rppa,countLambda=250,k=kvec[j])
  edgeset[[j]] = network.test.edges(ridgefit$pcor,verbose=FALSE,plot=FALSE)
}#end for j

# Remove list elements that have 0 rows (no edges predicted)
names(edgeset) <- paste0("k",kvec)
keep.ix <- which(lapply(edgeset,nrow) > 0)
edgeset <- edgeset[keep.ix]
saveRDS(edgeset,paste0(homedir,"/OUTPUT/EDGELISTS",network.method,".rds"))

#################################################
### LASSONET
#################################################
network.method = "lassonet"
# Vary parameter lambda
# 0 not included, causes convergence problems
lambdaVec = setdiff(seq(0,1,length.out=150),0)
edgeset = vector("list",length(lambdaVec))

for(j in 1:length(lambdaVec)){
  print(lambdaVec[j])
  # Get matrix of weights
  res <- lasso.net.vLambda(rppa,lambda=lambdaVec[j])
  print(res$lambda)
  edges = network.test.edges(res$pcor,verbose=FALSE,plot=FALSE)
  edgeset[[j]] = edges[which(abs(edges[,1]) > 0),]
}#end for j

# Remove list elements that have 0 rows (no edges predicted)
names(edgeset) <- paste0("lambda",lambdaVec)
keep.ix <- which(lapply(edgeset,nrow) > 0)
edgeset <- edgeset[keep.ix]
saveRDS(edgeset,paste0(homedir,"/OUTPUT/EDGELISTS",network.method,".rds"))

#################################################
### ELASTICNET
### The elasticNetwork method optimizes alpha and k internally using a grid search
#################################################
alphaVec <- seq(0.01,0.99,length.out<-10)
kvec <- c(3,4,5,10,20)

elfit <- elasticNetwork(rppa,nfold<-kvec,alpha<-alphaVec,verbose<-FALSE)
edges <- network.test.edges(elfit$pcor,verbose<-FALSE,plot<-FALSE)
edges$node1name <- prots[edges$node1]
edges$node2name <- prots[edges$node2]
write.table(edges,file.path(homedir,"OUTPUT/predictions_elasticnet.csv"),sep<-",",quote<-FALSE,row.names<-FALSE)

##################################################################################################
# Mutual information based methods (CLR, MRNET, ARACNE_additive, ARACNE_multiplicative)
##################################################################################################

#################################################
### CLR and MRNET
#################################################
network.method  <-  c("CLR","MRNET")

# Vary k in k-nearest neighbors
knnvec <- c(2,3,4,5,6)

edgesetCLR <- edgesetMRNET <- vector("list",length(knnvec))

for(j in 1:length(knnvec)){
  print(knnvec[j])
  mi <- knnmi.all(t(rppa),k = knnvec[j])
  
  ### CLR
  res <- clr(mi)
  edges <- getEdges(res)
  edgesetCLR[[j]] <- edges[,c(3,1,2)]
  
  ### MRNET
  res <- mrnet(mi)
  edges <- getEdges(res)
  edgesetMRNET[[j]] <- edges[,c(3,1,2)]
}#end for j

names(edgesetMRNET) <- names(edgesetCLR) <- paste0("k",knnvec)

# Remove list elements that have 0 rows (no edges predicted)
keep.ix <- which(lapply(edgesetCLR,nrow) > 0)
edgesetCLR <- edgesetCLR[keep.ix]
saveRDS(edgesetCLR,paste0(homedir,"/OUTPUT/EDGELISTS/",network.method[1],".rds"))  

# Remove list elements that have 0 rows (no edges predicted)
keep.ix <- which(lapply(edgesetMRNET,nrow) > 0)
edgesetMRNET <- edgesetMRNET[keep.ix]
saveRDS(edgesetMRNET,paste0(homedir,"/OUTPUT/EDGELISTS/",network.method[2],".rds"))  

#################################################
### ARACNE_A and ARACNE_M
#################################################
network.method <- c("aracne_a","aracne_m")

# Vary k in k-nearest neighbors
knnvec <- c(2,3,4,5,6)

# Several small values are added to epsilon and tau vectors because they perform best in some conditions
jnk1 <- c(1e-5,1e-4,seq(0,0.2,length.out<-100),seq(0.2,0.4,length.out = 49))
epsvec <- jnk1[-c(3,102)] # remove 0 and duplicates
jnk2 <- c(1e-5,1e-4,seq(0,0.5,length.out<-100),seq(0.5,1,length.out = 49))
tauvec <- jnk2[-c(3,102)] # remove 0 and duplicates

edgesetA <- edgesetM <- vector("list",length(knnvec))
names(edgesetA) <- names(edgesetM) <- paste0("k",knnvec)

for(j in 1:length(knnvec)){
  print(knnvec[j])
  mi <- knnmi.all(t(rppa),k=knnvec[j])
  edgesetA[[j]] <- edgesetM[[j]] <- vector("list",length(epsvec))
  names(edgesetA[[j]]) <- paste0("eps",epsvec)
  names(edgesetM[[j]]) <- paste0("tau",tauvec)
  for(z in 1:length(epsvec)){ 
    print(z)
    ### ARACNE_A
    res <- aracne.a(mi,eps = epsvec[z])
    edges <- getEdges(res)
    edgesetA[[j]][[z]] <- edges[,c(3,1,2)]
    
    ### ARACNE_M
    res <- aracne.m(mi,tau = tauvec[z])
    edges <- getEdges(res)
    edgesetM[[j]][[z]] <- edges[,c(3,1,2)]
    
  }#end for z
}#end for j

saveRDS(edgesetA,paste0(homedir,"/OUTPUT/EDGELISTS/",network.method[1],".rds"))
saveRDS(edgesetM,paste0(homedir,"/OUTPUT/EDGELISTS/",network.method[2],".rds"))  


