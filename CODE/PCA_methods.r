##################################################################################################
# Principal component analysis (PCA) on 13 network inference methods
##################################################################################################
# Dependencies for variables
# MAIN.r
##################################################################################################
methods = c("genenet","spearmancor", "pearsoncor", "simpleparcor","glasso","plsnet",
            "ridgenet","lassonet","elasticnet","aracne_a","aracne_m","clr_normalized","mrnet")

M = c()
for(i in 1:length(types)){
  print(types[i])
  # Find the union of genes for each tumor type
  unionset = c()
  for(j in 1:length(methods)){
    x = read.csv(paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_",methods[j],".txt",sep=""),sep="\t",header=T)
    edges = paste(x[,2],x[,3],sep="&")
    unionset = union(unionset,edges)
  }#end for j
  
  # Create a matrix to hold weights
  mat = matrix(0,nrow=length(unionset),ncol=length(methods))
  rownames(mat) = unionset
  colnames(mat) = methods
  for(j in 1:length(methods)){
    x = read.csv(paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_",methods[j],".txt",sep=""),sep="\t",header=T)
    edges = paste(x[,2],x[,3],sep="&")
    mx = match(edges,rownames(mat))
    mat[mx,j] = x[,1]
  }#end for j
  
  tmp = cbind(TUMOR=rep(types[i],nrow(mat)),mat)
  M = rbind(M,tmp)
}#end for i
write.table(M,"matrix_M_60202rows_13methods.txt",sep="\t",quote=F)

# PCA on the METHOD columns (remove the TUMOR column) yields the graphs in Fig 3d.

##################################################################################################
##################################################################################################