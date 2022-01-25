#################################################
# Run 13 network inference methods
#################################################
# Dependencies for variables
# MAIN.r
#################################################
# Required libraries and functions
#################################################
library(parcor)
library(glasso)
library(parmigene) # mutual information methods
source("functions.r")
#################################################
# GeneNet
#################################################
time = rep(NA,length(types))
### Read in RPPA files and get genenet p-values for all edges
for(i in 1:length(types)){
  print(i)
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  
  # GeneNet
  time[i] = system.time(res <- ggm.estimate.pcor(rppa))[3]
  all.edges = ggm.test.edges(res,verbose=FALSE,plot=FALSE)
  saveRDS(all.edges,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_genenet_from_pancan.rdata",sep=""))  
}#end for i
saveRDS(time,file="DATA/RUNNING TIME/genenet_time_all_tumors.rdata")

#################################################
### Spearman correlation
#################################################

time = rep(NA,length(types))
### Read in RPPA files and get genenet p-values for all edges
for(i in 1:length(types)){
  print(i)
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  # Get matrix of weights
  time[i] = system.time(res <- cor(rppa,method="spearman"))[3]
  all.edges = ggm.test.edges(res,verbose=FALSE,plot=FALSE)
  saveRDS(all.edges,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_spearmanCor_from_pancan.rdata",sep=""))  
}#end for i
saveRDS(time,file="DATA/RUNNING TIME/spearmanCor_time_all_tumors.rdata")

#################################################
# Pearson correlation
#################################################

time = rep(NA,length(types))
### Read in RPPA files and get genenet p-values for all edges
for(i in 1:length(types)){
  print(i)
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  # Get matrix of weights
  time[i] = system.time(res <- cor(rppa,method="pearson"))[3]
  all.edges = ggm.test.edges(res,verbose=FALSE,plot=FALSE)
  saveRDS(all.edges,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_pearsonCor_from_pancan.rdata",sep=""))  
}#end for i
saveRDS(time,file="DATA/RUNNING TIME/pearsonCor_time_all_tumors.rdata")

#################################################
# Simple partial correlation
#################################################

time = rep(NA,length(types))
### Read in RPPA files and get genenet p-values for all edges
for(i in 1:length(types)){
  print(i)
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  # Get matrix of weights
  time[i] = system.time(res <- cor2pcor(cor(rppa,method="pearson")))[3]
  all.edges = ggm.test.edges(res,verbose=FALSE,plot=FALSE)
  saveRDS(all.edges,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_simpleParCor_from_pancan.rdata",sep=""))  
}#end for i
saveRDS(time,file="DATA/RUNNING TIME/simpleParCor_time_all_tumors.rdata")

#################################################
### GLASSO
#################################################
# GLASSO auxiliary function
runGlasso = function(covrppa,rho=0.01){
  gl <- glasso(covrppa,rho=rho)
  res <- invcov2parcor(gl$wi)
  return(res)
}#end function

jnk = seq(0,1,length.out=150)
# 0 causes convergence problems, so we take that out
rhoVec = jnk[-1] # 149 entries

times = matrix(NA,nrow=length(rhoVec),ncol=length(types))
colnames(times) = types
### Read in RPPA files and get p-values for all edges
for(i in 1:length(types)){
  print(i)
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  edgeset = vector("list",length(rhoVec))
  covrppa = cov(rppa)
  
  for(j in 1:length(rhoVec)){
    print(j)
    # Get matrix of weights
    times[j,i] = system.time(res <- runGlasso(covrppa,rho=rhoVec[j]))[3]
    all.edges = ggm.test.edges(res,verbose=FALSE,plot=FALSE)
    edgeset[[j]] = all.edges[which(abs(all.edges[,1]) > 0),]
  }#end for j
  saveRDS(edgeset,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_glasso_from_pancan.rdata",sep=""))  
}#end for i
saveRDS(times,file="DATA/RUNNING TIME/glasso_time_all_tumors.rdata")

#################################################
### Partial Least Squares
#################################################

kvec = c(3,4,5,10,20)
times = matrix(NA,nrow=length(kvec),ncol=length(types))
colnames(times) = types
rownames(times) = paste("K",kvec,sep="")
ncomp = 30
network.method = "plsnet"

for(i in 1:length(types)){
  print(i)
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  rppa = as.matrix(rppa)
  edgeset = vector("list",length(kvec))
  for(j in 1:length(kvec)){
    cat("j is",j,"\n")
    times[j,i] = system.time(plsfit <- pls.net(rppa,scale=TRUE,k=kvec[j],ncomp=30,verbose=FALSE))[3]
    edgeset[[j]] = ggm.test.edges(plsfit$pcor,verbose=FALSE,plot=FALSE)
    print(plsfit$m)
  }#end for j
  saveRDS(edgeset,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_plsnet_from_pancan.rdata",sep=""))  
}#end for i
saveRDS(times,file="DATA/RUNNING TIME/plsnet_time_all_tumors.rdata")

#################################################
### RIDGENET
#################################################

kvec = c(3,4,5,10,20)
times = matrix(NA,nrow=length(kvec),ncol=length(types))
colnames(times) = types
rownames(times) = paste("K",kvec,sep="")
network.method = "ridgenet"

for(i in 1:length(types)){
  print(i)
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  rppa = as.matrix(rppa)
  edgeset = vector("list",length(kvec))
  for(j in 1:length(kvec)){
    cat("j is",j,"\n")
    times[j,i] = system.time(ridgefit <- ridge.net.vLambda(rppa,countLambda=250,k=kvec[j]))[3]
    edgeset[[j]] = ggm.test.edges(ridgefit$pcor,verbose=FALSE,plot=FALSE)
  }#end for j
  saveRDS(edgeset,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_ridgenet_from_pancan.rdata",sep=""))  
}#end for i
saveRDS(times,file="DATA/RUNNING TIME/ridgenet_time_all_tumors.rdata")

#################################################
### LASSONET
#################################################

network.method = "lassonet"

jnk = seq(0,1,length.out=150)
# 0 causes convergence problems, so we take that out
lambdaVec = jnk[-1] # 149 entries

times = matrix(NA,nrow=length(lambdaVec),ncol=length(types))
colnames(times) = types
### Read in RPPA files and get p-values for all edges
for(i in 1:length(types)){
  print(i)
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  edgeset = vector("list",length(lambdaVec))
  
  for(j in 1:length(lambdaVec)){
    print(lambdaVec[j])
    # Get matrix of weights
    times[j,i] = system.time(res <- lasso.net.vLambda(rppa,lambda=lambdaVec[j]))[3]
    print(res$lambda)
    all.edges = ggm.test.edges(res$pcor,verbose=FALSE,plot=FALSE)
    edgeset[[j]] = all.edges[which(abs(all.edges[,1]) > 0),]
  }#end for j
  saveRDS(edgeset,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_lassonet_from_pancan.rdata",sep=""))  
}#end for i
saveRDS(times,file="DATA/RUNNING TIME/lassonet_time_all_tumors.rdata")

#################################################
### ELASTICNET
#################################################
# alpha will be optimized together with k using grid search

alphaVec = seq(0.01,0.99,length.out=10)
kvec = c(3,4,5,10,20)

times = rep(NA,length(types))
names(times) = types
network.method = "elasticnet"
opt.alpha = opt.K = matrix(NA,nrow=ncol(rppa),ncol=length(types))
colnames(opt.alpha) = colnames(opt.K) = types
rownames(opt.alpha) = rownames(opt.K) = 1:ncol(rppa)

for(i in 1:length(types)){
  print(i)
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  rppa = as.matrix(rppa)
  times[i] = system.time(elfit <- elasticNetwork(rppa,nfold=kvec,alpha=alphaVec,verbose=FALSE))[3]
  edgeset = ggm.test.edges(elfit$pcor,verbose=FALSE,plot=FALSE)
  opt.alpha[,i] = elfit$alpha.opt
  opt.K[,i] = elfit$n.opt
  saveRDS(edgeset,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_elasticnet_from_pancan.rdata",sep="")) 
}#end for i
saveRDS(times,file="DATA/RUNNING TIME/elasticnet_time_all_tumors.rdata")

#################################################
### CLR and MRNET
#################################################

knnvec = c(2,3,4,5,6)
timesCLR = timesMRNET = matrix(NA,nrow=length(knnvec),ncol=length(types)) 
rownames(timesCLR) = rownames(timesMRNET) = paste("knn",knnvec,sep="")
colnames(timesCLR) = colnames(timesMRNET) = types

for(i in 1:length(types)){
  print(types[i])
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  
  edgesetCLR = edgesetMRNET = vector("list",length(knnvec))
  
  for(j in 1:length(knnvec)){
    print(knnvec[j])
    mitime = system.time(mi <- knnmi.all(t(rppa),k=knnvec[j]))[3]
    
    ### CLR
    timesCLR[j,i] = system.time(res <- clr(mi))[3] + mitime
    all.edges = getEdges(res)
    edgesetCLR[[j]] = all.edges[,c(3,1,2)]
    
    ### MRNET
    timesMRNET[j,i] = system.time(res <- mrnet(mi))[3] + mitime
    all.edges = getEdges(res)
    edgesetMRNET[[j]] = all.edges[,c(3,1,2)]
  }#end for j
  
  saveRDS(edgesetCLR,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_clr_from_pancan.rdata",sep=""))  
  saveRDS(edgesetMRNET,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_mrnet_from_pancan.rdata",sep=""))  
  
}#end for i

saveRDS(timesCLR,file=paste("DATA/RUNNING TIME/clr_time_all_tumors.rdata",sep=""))
saveRDS(timesMRNET,file=paste("DATA/RUNNING TIME/mrnet_time_all_tumors.rdata",sep=""))

#################################################
### ARACNE_A and ARACNE_M
#################################################

network.method = "aracne"
knnvec = c(2,3,4,5,6)
jnk1 = c(1e-5,1e-4,seq(0,0.2,length.out=100),seq(0.2,0.4,length.out=49))
epsvec = jnk1[-c(3,102)]
jnk2 = c(1e-5,1e-4,seq(0,0.5,length.out=100),seq(0.5,1,length.out=49))
tauvec = jnk2[-c(3,102)]


### Read in RPPA files and get p-values for all edges
for(i in 1:length(types)){
  print(types[i])
  rppa = readRDS(paste("DATA/INPUT RPPA DATA/useable_",types[i],"_from_pancan.rdata",sep=""))
  
  mitime = rep(NA,length(knnvec)) # running time is the sum of mutual information time and aracne time 
  timesA = timesM = matrix(NA,nrow=length(epsvec),ncol=length(knnvec))
  colnames(timesA) = colnames(timesM) = paste("knn",knnvec,sep="")
  
  for(j in 1:length(knnvec)){
    print(knnvec[j])
    mitime[j] = system.time(mi <- knnmi.all(t(rppa),k=knnvec[j]))[3]
    edgesetA = edgesetM = vector("list",length(epsvec))
    
    for(z in 1:length(epsvec)){ 
      print(z)
      ### ARACNE_A
      timesA[z,j] = system.time(res <- aracne.a(mi,eps=epsvec[z]))[3]
      all.edges = getEdges(res)
      edgesetA[[z]] = all.edges[,c(3,1,2)]
      
      ### ARACNE_M
      timesM[z,j] = system.time(res <- aracne.m(mi,tau=tauvec[z]))[3]
      all.edges = getEdges(res)
      edgesetM[[z]] = all.edges[,c(3,1,2)]
      
    }#end for z
    saveRDS(edgesetA,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_aracne_a_knn",knnvec[j],"from_pancan.rdata",sep=""))  
    saveRDS(edgesetM,file=paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_aracne_m_knn",knnvec[j],"from_pancan.rdata",sep=""))  
    
    timesA[,j] = timesA[,j] + mitime[j]
    timesM[,j] = timesM[,j] + mitime[j]
  }#end for j
  saveRDS(timesA,file=paste("DATA/RUNNING TIME/aracne_a_time_",types[i],".rdata",sep=""))
  saveRDS(timesM,file=paste("DATA/RUNNING TIME/aracne_m_time_",types[i],".rdata",sep=""))
}#end for i

##################################################################################################
##################################################################################################

