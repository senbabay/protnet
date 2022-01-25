##################################################################################################
# Normalizing CLR weights
##################################################################################################
# Dependencies for variables
# MAIN.r
##################################################################################################
a = readRDS("DATA/LIMITED RECALL OPTIMAL PARAMETER VALUES/clr_AUPR_limitedRecall_best_K.rdata")
for(i in 1:length(types)){
  edgelist = readRDS(paste("DATA/OUTPUT EDGELISTS/edges_",types[i],"_clr_from_pancan.rdata",sep=""))
  j = a[i,1]
  x = edgelist[[j]]
  normalized = x[,1]/sum(x[,1])
  NORM_CONSTANT = 1 / max(normalized) # normalization constant that brings the range of coefficients to [0,1]
  newCoefs = normalized * NORM_CONSTANT
  dot1 = apply(as.matrix(x[,2:3]),1,paste,collapse=".")
  
  edgefile = read.csv(paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_clr.txt",sep=""),sep="\t",header=T)
  dot2 = apply(as.matrix(edgefile[,2:3]),1,paste,collapse=".")
  mx = match(dot2,dot1)
  newfile = edgefile
  newfile[,1] = newCoefs[mx]
  
  write.table(newfile,paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_clr_normalized.txt",sep=""),sep="\t",quote=F,row.names=F)
}#end