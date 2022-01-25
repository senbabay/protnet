##################################################################################################
# Choosing the best parameter value for the 9 methods that require user optimization
# Results are reported in S5 Fig.
##################################################################################################
# Dependencies for variables
# MAIN.r
##################################################################################################
# GLASSO
##################################################################################################
areaVec = matrix(NA,nrow=length(types),ncol=149)
sortedEdges = vector("list",149)
bestEdges = vector("list",length(types))

for(i in 1:length(types)){
  print(types[i])
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_glasso_from_pancan.rdata",sep=""))
  for(j in 1:length(edgelist)){
    if(j %% 20 == 0){
      print(j)
    }
    res = longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE)
    areaVec[i,j] = res$area
    sortedEdges[[j]] = res$sortedEdges
  }#end for j
  bestEdges[[i]] = sortedEdges[[which.max(areaVec[i,])]]
}#end for i

jnk = seq(0,1,length.out=150)
rhoVec = jnk[-1] # 149 entries

areaVal = apply(areaVec,1,max)
bestRhoInd = apply(areaVec,1,which.max)
bestRho = rhoVec[bestRhoInd]
newArea = matrix(NA,nrow=length(types),ncol=3)
rownames(newArea) = types
colnames(newArea) = c("best.rho.index","best.rho.val","best.rho.AUPR")
newArea[,1] = bestRhoInd
newArea[,2] = bestRho
newArea[,3] = areaVal

saveRDS(newArea,"DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/glasso_AUPR_limitedRecall_best_rho.rdata")
saveRDS(bestEdges,"DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/glasso_EDGES_limitedRecall_best_rho.rdata")
##################################################################################################
# PLSNET
##################################################################################################

areaVec = matrix(NA,nrow=length(types),ncol=5)
kvec = c(3,4,5,10,20)
sortedEdges = vector("list",5)
bestEdges = vector("list",length(types))

for(i in 1:length(types)){
  print(types[i])
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_plsnet_from_pancan.rdata",sep=""))
  for(j in 1:length(edgelist)){
    res = longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE)
    areaVec[i,j] = res$area
    sortedEdges[[j]] = res$sortedEdges
  }#end for j
  bestEdges[[i]] = sortedEdges[[which.max(areaVec[i,])]]
}#end for i

areaVal = apply(areaVec,1,max)
bestKInd = apply(areaVec,1,which.max)
bestK = kvec[bestKInd]

newArea = matrix(NA,nrow=length(types),ncol=3)
rownames(newArea) = types
colnames(newArea) = c("best.K.index","best.K.val","best.K.AUPR")
newArea[,1] = bestKInd
newArea[,2] = bestK
newArea[,3] = areaVal

saveRDS(newArea,"DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/plsnet_AUPR_limitedRecall_best_K.rdata")
saveRDS(bestEdges,"DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/plsnet_EDGES_limitedRecall_best_K.rdata")

##################################################################################################
# RIDGENET
##################################################################################################
areaVec = matrix(NA,nrow=length(types),ncol=5)
kvec = c(3,4,5,10,20)
sortedEdges = vector("list",5)
bestEdges = vector("list",length(types))

for(i in 1:length(types)){
  print(types[i])
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_ridgenet_from_pancan.rdata",sep=""))
  for(j in 1:length(edgelist)){
    res = longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE)
    areaVec[i,j] = res$area
    sortedEdges[[j]] = res$sortedEdges
  }#end for j
}#end for i

areaVal = apply(areaVec,1,max)
bestKInd = apply(areaVec,1,which.max)
bestK = kvec[bestKInd]

newArea = matrix(NA,nrow=length(types),ncol=3)
rownames(newArea) = types
colnames(newArea) = c("best.K.index","best.K.val","best.K.AUPR")
newArea[,1] = bestKInd
newArea[,2] = bestK
newArea[,3] = areaVal

saveRDS(newArea,"DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/ridgenet_AUPR_limitedRecall_best_K.rdata")
saveRDS(bestEdges,"DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/ridgenet_EDGES_limitedRecall_best_K.rdata")

##################################################################################################
# LASSONET
##################################################################################################

areaVec = matrix(NA,nrow=length(types),ncol=149)
sortedEdges = vector("list",149)
bestEdges = vector("list",length(types))

jnk = seq(0,1,length.out=150)
lambdaVec = jnk[-1] # 149 entries

for(i in 1:length(types)){
  print(types[i])
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_lassonet_from_pancan.rdata",sep=""))
  for(j in 1:length(edgelist)){
    if(j %% 20 == 0){
      print(j)
    }
    res = longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE)
    areaVec[i,j] = res$area
    sortedEdges[[j]] = res$sortedEdges
  }#end for j
  bestEdges[[i]] = sortedEdges[[which.max(areaVec[i,])]]
}#end for i

areaVal = apply(areaVec,1,max)
bestLInd = apply(areaVec,1,which.max)
bestL = lambdaVec[bestLInd]
newArea = matrix(NA,nrow=length(types),ncol=3)
rownames(newArea) = types
colnames(newArea) = c("best.lambda.index","best.lambda.val","best.lambda.AUPR")
newArea[,1] = bestLInd
newArea[,2] = bestL
newArea[,3] = areaVal

saveRDS(newArea,"DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/lassonet_AUPR_limitedRecall_best_lambda.rdata")
saveRDS(bestEdges,"DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/lassonet_EDGES_limitedRecall_best_lambda.rdata")

##################################################################################################
# ELASTICNET
##################################################################################################

areaVec = matrix(NA,nrow=length(types),ncol=1)
colnames(areaVec) = "AUPR"
rownames(areaVec) = types
bestEdges = vector("list",length(types))

for(i in 1:length(types)){
  print(types[i])
  edgefile = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_elasticnet_from_pancan.rdata",sep=""))
  last2keep = min(which(edgefile[[1]] == 0)) - 1
  edgelist = cbind(edgefile[[2]][1:last2keep],edgefile[[3]][1:last2keep])
  res = longPrecisionRecall(edgelist,pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE,node1ind=1,node2ind=2)
  areaVec[i] = res$area
  bestEdges[[i]] = res$sortedEdges
}#end for i

saveRDS(areaVec,"DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/elasticnet_AUPR_limitedRecall_best_alphaAndK.rdata")
saveRDS(bestEdges,"DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/elasticnet_EDGES_limitedRecall_best_alphaAndK.rdata")

##################################################################################################
# ARACNE.A
##################################################################################################
knnvec = c(2,3,4,5,6)
jnk1 = c(1e-5,1e-4,seq(0,0.2,length.out=100),seq(0.2,0.4,length.out=49))
epsvec = jnk1[-c(3,102)]
# Optimize k and epsilon together (grid search)
# meaning find the maximum area in the matrix of k and epsilon
bestK = bestEps = bestArea = rep(NA,length(types))
bestEdges = vector("list",length(types))

for(i in 1:length(types)){
  print(types[i])
  maxArea = 0
  for(j in 1:length(knnvec)){
    print(knnvec[j])
    edgefile = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_aracne_a_knn",knnvec[j],"from_pancan.rdata",sep=""))
    for(w in 1:length(edgefile)){
      if(w %% 20 == 0){
        print(w)
      }
      res = longPrecisionRecall(edgefile[[w]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE,node1ind=2,node2ind=3)
      if(res$area > maxArea){
        maxArea = res$area
        bestArea[i] = res$area
        bestK[i] = knnvec[j]
        bestEps[i] = epsvec[w]
        bestEdges[[i]] = res$sortedEdges
      }#end if
    }#end for w
  }#end for j
}#end for i

areaVec = cbind(best.K.val=bestK,best.eps.val=bestEps,best.AUPR=bestArea)
rownames(areaVec) = types
saveRDS(areaVec,"DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/aracne_a_AUPR_limitedRecall_best_EpsAndK.rdata")
saveRDS(bestEdges,"DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/aracne_a_EDGES_limitedRecall_best_EpsAndK.rdata")

##################################################################################################
# ARACNE.M
##################################################################################################
jnk2 = c(1e-5,1e-4,seq(0,0.5,length.out=100),seq(0.5,1,length.out=49))
tauvec = jnk2[-c(3,102)]
# Optimize k and tau together (grid search)
# meaning find the maximum area in the matrix of tau and epsilon
bestK = bestTau = bestArea = rep(NA,length(types))
bestEdges = vector("list",length(types))

for(i in 1:length(types)){
  print(types[i])
  maxArea = 0
  for(j in 1:length(knnvec)){
    print(knnvec[j])
    edgefile = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_aracne_m_knn",knnvec[j],"from_pancan.rdata",sep=""))
    for(w in 1:length(edgefile)){
      if(w %% 20 == 0){
        print(w)
      }
      res = longPrecisionRecall(edgefile[[w]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE,node1ind=2,node2ind=3)
      if(res$area > maxArea){
        maxArea = res$area
        bestArea[i] = res$area
        bestK[i] = knnvec[j]
        bestTau[i] = tauvec[w]
        bestEdges[[i]] = res$sortedEdges
      }#end if
    }#end for w
  }#end for j
}#end for i

areaVec = cbind(best.K.val=bestK,best.tau.val=bestTau,best.AUPR=bestArea)
rownames(areaVec) = types
saveRDS(areaVec,"DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/aracne_m_AUPR_limitedRecall_best_TauAndK.rdata")
saveRDS(bestEdges,"DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/aracne_m_EDGES_limitedRecall_best_TauAndK.rdata")

##################################################################################################
# CLR
##################################################################################################
areaVec = matrix(NA,nrow=length(types),ncol=5)
knnvec = c(2,3,4,5,6)
sortedEdges = vector("list",5)
bestEdges = vector("list",length(types))

for(i in 1:length(types)){
  print(types[i])
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_clr_from_pancan.rdata",sep=""))
  for(j in 1:length(edgelist)){
    res = longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE)
    areaVec[i,j] = res$area
    sortedEdges[[j]] = res$sortedEdges
  }#end for j
  bestEdges[[i]] = sortedEdges[[which.max(areaVec[i,])]]
}#end for i

areaVal = apply(areaVec,1,max)
bestKInd = apply(areaVec,1,which.max)
bestK = knnvec[bestKInd]

newArea = matrix(NA,nrow=length(types),ncol=3)
rownames(newArea) = types
colnames(newArea) = c("best.K.index","best.K.val","best.K.AUPR")
newArea[,1] = bestKInd
newArea[,2] = bestK
newArea[,3] = areaVal

saveRDS(newArea,"DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/clr_AUPR_limitedRecall_best_K.rdata")
saveRDS(bestEdges,"DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/clr_EDGES_limitedRecall_best_K.rdata")

##################################################################################################
# MRNET
##################################################################################################
areaVec = matrix(NA,nrow=length(types),ncol=5)
knnvec = c(2,3,4,5,6)
sortedEdges = vector("list",5)
bestEdges = vector("list",length(types))

for(i in 1:length(types)){
  print(types[i])
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_mrnet_from_pancan.rdata",sep=""))
  for(j in 1:length(edgelist)){
    res = longPrecisionRecall(edgelist[[j]],pcInts,prots=prots,recall.limit=0.1,sparsePred=TRUE)
    areaVec[i,j] = res$area
    sortedEdges[[j]] = res$sortedEdges
  }#end for j
  bestEdges[[i]] = sortedEdges[[which.max(areaVec[i,])]]
  #plot(res$recall,res$prec,type="l",col=4)
}#end for i

areaVal = apply(areaVec,1,max)
bestKInd = apply(areaVec,1,which.max)
bestK = knnvec[bestKInd]

newArea = matrix(NA,nrow=length(types),ncol=3)
rownames(newArea) = types
colnames(newArea) = c("best.K.index","best.K.val","best.K.AUPR")
newArea[,1] = bestKInd
newArea[,2] = bestK
newArea[,3] = areaVal

saveRDS(newArea,paste("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/mrnet_AUPR_limitedRecall_best_K.rdata",sep=""))
saveRDS(bestEdges,paste("DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/mrnet_EDGES_limitedRecall_best_K.rdata",sep=""))
##################################################################################################
##################################################################################################







