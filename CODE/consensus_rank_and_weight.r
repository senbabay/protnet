##################################################################################################
# Consensus edge ranks and edge weights across top 6 methods
##################################################################################################
# Dependencies for variables
# MAIN.r, PCA_methods.r
##################################################################################################

allPairs = t(combn(1:187,2))
allPossibleEdges = paste(allPairs[,1],allPairs[,2],sep="&") # N=17391

# Top 6 methods "spearmancor"    "ridgenet"       "lassonet"       "aracne_a"       "aracne_m"       "clr"
ix = c(2,7,8,10,11,12)
methods.top6 = methods[ix]

for(i in 1:length(types)){
  print(types[i])
  # Create a master matrix of ranks with all possible edges, and then insert the ranks from each method
  ranks = matrix(NA,nrow=length(allPossibleEdges),ncol=length(methods.top6)) 
  rownames(ranks) = allPossibleEdges
  colnames(ranks) = methods.top6
  # Also create a master matrix of weights
  weights = ranks
  
  for(j in 1:length(methods.top6)){
    # Limited recall edges
    x = read.csv(paste("DATA/LIMITED RECALL EDGES AND WEIGHTS/edges_and_weights_limited_recall_",types[i],"_",methods.top6[j],".txt",sep=""),sep="\t",header=T)
    # this.rank = rank(-x[,1],ties.method = "min")
    
    # Note that edges that are not in X get rank [1+NUMBER OF NON-ZERO EDGES IN X]
    # We do not include the zero edges in ranking, because they were randomly chosen
    # so their ranks cannot be used as a proxy for significance
    
    # So make [1+NUMBER OF NON-ZERO EDGES IN X] the default rank for all edges in the relevant column of RANKS
    x.cut = x[which(abs(x[,1]) > 0),]
    ranks[,j] = 1 + nrow(x.cut)
    
    ordx = order(as.numeric(abs(x.cut[,1])),decreasing=TRUE)
    x2 = x.cut[ordx,]
    
    # connect edges
    edges = paste(x2[,2],x2[,3],sep="&")
    
    # Find edges in the universe set
    mx = match(edges,allPossibleEdges)
    
    # Place inside RANKS the ranks of edges in X  
    ranks[mx,j] = 1:length(edges)
    
    # Place inside WEIGHTS the weights of edges in X
    # Keep NA weights for edges not in X
    weights[mx,j] = x2[,1]
  }#end for j
  
  # Average across top 6 methods
  aveRank = apply(ranks,1,mean) # does not have any NAs
  aveWeight = apply(abs(weights),1,mean,na.rm=TRUE) # remove NAs
  
  # Remove entries where average weight is NaN (all methods had NA for this edge) or 0 (from one of the sparse methods)
  w = which(aveWeight=="NaN" | aveWeight==0)
  aveRank2 = aveRank[-w]
  aveWeight2 = aveWeight[-w]
  rm(aveRank,aveWeight)
  
  # Sign of the edges will be derived from the sign of the SPEARMAN correlations
  protMat = pancan[which(pancan[,2]==types[i]),-c(1,2)]
  signs = sign(cor(protMat,method="spearman"))
  signvec = signs[lower.tri(signs)]
  names(signvec) = allPossibleEdges
  mx = match(names(aveWeight2),allPossibleEdges)
  spearmanSign = signvec[mx]
  
  outmat = cbind(EDGE=names(aveWeight2),consensusRank=round(as.numeric(aveRank2),3),consensusWeight=round(as.numeric(aveWeight2),5),spearmanSign=spearmanSign)
  # Order according to consensusRank
  op = order(aveRank2)
  outmat2 = outmat[op,]
  rm(outmat)
  
  write.table(outmat2,file=paste("DATA/CONSENSUS RANKS AND WEIGHTS/top6Methods_",types[i],".txt",sep=""),quote=F,sep="\t",row.names=F)
  
}#end for i

