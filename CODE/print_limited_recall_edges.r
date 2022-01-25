##################################################################################################
# Printing edges and their weights from the "limited recall" run
##################################################################################################
# Dependencies for variables
# MAIN.r, choose_best_param.r
##################################################################################################
# GLASSO
##################################################################################################
# Get number of edges at limited recall, and create a weight matrix including only these edges
a = readRDS("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/glasso_AUPR_limitedRecall_best_rho.rdata")
b = readRDS("DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/glasso_EDGES_limitedRecall_best_rho.rdata")
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_glasso_from_pancan.rdata",sep=""))
  m = longPrecisionRecall(b[[i]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[[i]][1:count,]
  weights = rep(0,nrow(edges))
  
  j = a[i,1]
  # edgelist[[j]] gives you weights for the top X, the remaining get weight 0
  x = edgelist[[j]]
  np = min(nrow(x),length(weights))
  weights[1:np] = -x[1:np,1] # GLASSO coefficients have to be multiplied by 1
  outmat = cbind(weights,edges)
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_glasso.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end

##################################################################################################
# PLSNET
##################################################################################################
# Get number of edges at limited recall, and create a weight matrix including only these edges
a = readRDS("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/plsnet_AUPR_limitedRecall_best_K.rdata")
b = readRDS("DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/plsnet_EDGES_limitedRecall_best_K.rdata")
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_plsnet_from_pancan.rdata",sep=""))
  m = longPrecisionRecall(b[[i]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[[i]][1:count,]
  weights = rep(0,nrow(edges))
  
  j = a[i,1]
  x = edgelist[[j]]
  np = min(nrow(x),length(weights))
  weights[1:np] = x[1:np,1] 
  outmat = cbind(weights,edges)
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  print(nrow(outmat))
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_plsnet.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end

##################################################################################################
# RIDGENET
##################################################################################################
# Get number of edges at limited recall, and create a weight matrix including only these edges
a = readRDS("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/ridgenet_AUPR_limitedRecall_best_K.rdata")
b = readRDS("DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/ridgenet_EDGES_limitedRecall_best_K.rdata")
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_ridgenet_from_pancan.rdata",sep=""))
  m = longPrecisionRecall(b[[i]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[[i]][1:count,]
  weights = rep(0,nrow(edges))
  
  j = a[i,1]
  x = edgelist[[j]]
  np = min(nrow(x),length(weights))
  weights[1:np] = x[1:np,1] 
  outmat = cbind(weights,edges)
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  print(nrow(outmat))
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_ridgenet.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end

##################################################################################################
# LASSONET
##################################################################################################
# Get number of edges at limited recall, and create a weight matrix including only these edges
a = readRDS("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/lassonet_AUPR_limitedRecall_best_lambda.rdata")
b = readRDS("DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/lassonet_EDGES_limitedRecall_best_lambda.rdata")
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_lassonet_from_pancan.rdata",sep=""))
  m = longPrecisionRecall(b[[i]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[[i]][1:count,]
  weights = rep(0,nrow(edges))
  
  j = a[i,1]
  x = edgelist[[j]]
  np = min(nrow(x),length(weights))
  weights[1:np] = x[1:np,1] 
  outmat = cbind(weights,edges)
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  print(nrow(outmat))
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_lassonet.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end

##################################################################################################
# ELASTICNET
##################################################################################################
# Get number of edges at limited recall, and create a weight matrix including only these edges
a = readRDS("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/elasticnet_AUPR_limitedRecall_best_alphaAndK.rdata")
b = readRDS("DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/elasticnet_EDGES_limitedRecall_best_alphaAndK.rdata")
for(i in 1:length(types)){
  m = longPrecisionRecall(b[[i]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[[i]][1:count,]
  weights = rep(0,nrow(edges))
  
  edgefile = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_elasticnet_from_pancan.rdata",sep=""))
  last2keep = min(which(edgefile[[1]] == 0)) - 1
  x = cbind(edgefile[[1]][1:last2keep],edgefile[[2]][1:last2keep],edgefile[[3]][1:last2keep])
  
  np = min(nrow(x),length(weights))
  weights[1:np] = x[1:np,1] 
  outmat = cbind(weights,edges)
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  print(nrow(outmat))
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_elasticnet.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end

##################################################################################################
# ARACNE.A
##################################################################################################
# Get number of edges at limited recall, and create a weight matrix including only these edges
a = readRDS("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/aracne_a_AUPR_limitedRecall_best_EpsAndK.rdata")
b = readRDS("DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/aracne_a_EDGES_limitedRecall_best_EpsAndK.rdata")
for(i in 1:length(types)){
  m = longPrecisionRecall(b[[i]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[[i]][1:count,]
  weights = rep(0,nrow(edges))
  
  j = a[i,1]
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_aracne_a_knn",j,"from_pancan.rdata",sep=""))
  w = which(epsvec==a[i,2])
  x = edgelist[[w]]
  np = min(nrow(x),length(weights))
  weights[1:np] = x[1:np,1] 
  outmat = cbind(weights,edges)
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  print(nrow(outmat))
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_aracne_a.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end

##################################################################################################
# ARACNE.M
##################################################################################################
# Get number of edges at limited recall, and create a weight matrix including only these edges
a = readRDS("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/aracne_m_AUPR_limitedRecall_best_TauAndK.rdata")
b = readRDS("DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/aracne_m_EDGES_limitedRecall_best_TauAndK.rdata")
for(i in 1:length(types)){
  m = longPrecisionRecall(b[[i]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[[i]][1:count,]
  weights = rep(0,nrow(edges))
  
  j = a[i,1]
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_aracne_m_knn",j,"from_pancan.rdata",sep=""))
  w = which(tauvec==a[i,2])
  x = edgelist[[w]]
  np = min(nrow(x),length(weights))
  weights[1:np] = x[1:np,1] 
  outmat = cbind(weights,edges)
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  print(nrow(outmat))
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_aracne_m.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end

##################################################################################################
# CLR
##################################################################################################
# Get number of edges at limited recall, and create a weight matrix including only these edges
a = readRDS("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/clr_AUPR_limitedRecall_best_K.rdata")
b = readRDS("DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/clr_EDGES_limitedRecall_best_K.rdata")
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_clr_from_pancan.rdata",sep=""))
  m = longPrecisionRecall(b[[i]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[[i]][1:count,]
  weights = rep(0,nrow(edges))
  
  j = a[i,1]
  x = edgelist[[j]]
  np = min(nrow(x),length(weights))
  weights[1:np] = x[1:np,1] 
  outmat = cbind(weights,edges)
  print(nrow(outmat))
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_clr.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end

##################################################################################################
# MRNET
##################################################################################################
# Get number of edges at limited recall, and create a weight matrix including only these edges
a = readRDS("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/mrnet_AUPR_limitedRecall_best_K.rdata")
b = readRDS("DATA/INPUT EDGELISTS FOR LONG_PRECISION_RECALL/mrnet_EDGES_limitedRecall_best_K.rdata")
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_mrnet_from_pancan.rdata",sep=""))
  m = longPrecisionRecall(b[[i]],pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[[i]][1:count,]
  weights = rep(0,nrow(edges))
  
  j = a[i,1]
  x = edgelist[[j]]
  np = min(nrow(x),length(weights))
  weights[1:np] = x[1:np,1] 
  outmat = cbind(weights,edges)
  print(nrow(outmat))
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_mrnet.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end

##################################################################################################
# GENENET
##################################################################################################
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_genenet_from_pancan.rdata",sep=""))
  b = cbind(edgelist[[2]],edgelist[[3]])
  m = longPrecisionRecall(b,pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[1:count,]
  weights = edgelist[[1]][1:count]
  
  outmat = cbind(weights,edges)
  print(nrow(outmat))
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_genenet.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end for i

##################################################################################################
# SPEARMANCOR
##################################################################################################
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_spearmancor_from_pancan.rdata",sep=""))
  b = cbind(edgelist[[2]],edgelist[[3]])
  m = longPrecisionRecall(b,pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[1:count,]
  weights = edgelist[[1]][1:count]
  
  outmat = cbind(weights,edges)
  print(nrow(outmat))
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_spearmancor.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end for i

##################################################################################################
# PEARSONCOR
##################################################################################################
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_pearsoncor_from_pancan.rdata",sep=""))
  b = cbind(edgelist[[2]],edgelist[[3]])
  m = longPrecisionRecall(b,pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[1:count,]
  weights = edgelist[[1]][1:count]
  
  outmat = cbind(weights,edges)
  print(nrow(outmat))
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_pearsoncor.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end for i

##################################################################################################
# SIMPLEPARCOR
##################################################################################################
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_simpleparcor_from_pancan.rdata",sep=""))
  b = cbind(edgelist[[2]],edgelist[[3]])
  m = longPrecisionRecall(b,pcInts,prots=prots,recall.limit=0.1,sparsePred=FALSE,node1ind = 1,node2ind = 2)
  count = length(m$recall)
  edges = b[1:count,]
  weights = edgelist[[1]][1:count]
  
  outmat = cbind(weights,edges)
  print(nrow(outmat))
  colnames(outmat) = c("WEIGHTS","ANTIBODY 1","ANTIBODY 2")
  write.table(outmat,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_simpleparcor.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end for i
##################################################################################################
##################################################################################################





