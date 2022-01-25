##################################################################################################
# Create matrix T, the discovery set including 1008 edges and 11 tumor types
##################################################################################################
# Dependencies for variables
# MAIN.r
##################################################################################################
aveRankCutoff = 425
unionSet = c()
indEdgeset = vector("list",length(types))
names(indEdgeset) = types

# Collect from 11 tumor types all edges that pass threshold
for(i in 1:length(types)){
  indat = read.csv(paste("DATA/CONSENSUS RANKS AND WEIGHTS/top6Methods_",types[i],".txt",sep=""),sep="\t",header=TRUE)
  ind = which(indat[,2] <= aveRankCutoff)  
  print(length(ind))
  thisEdgeset = as.character(indat[ind,1])
  unionSet = union(thisEdgeset,unionSet)
  
  # Also store individually in list 
  indEdgeset[[i]] = thisEdgeset
}#end for i
# Union comes to 1008 edges
######################################################################
# Get consensus weights and binary calls for the edges
######################################################################
weightMat = matrix(NA,nrow=length(unionSet),ncol=length(types))
colnames(weightMat) = types
rownames(weightMat) = unionSet
signifBinMat = weightMat

for(i in 1:length(types)){
  indat = read.csv(paste("DATA/CONSENSUS RANKS AND WEIGHTS/top6Methods_",types[i],".txt",sep=""),sep="\t",header=TRUE,colClasses=c("character","numeric","numeric","integer"))
  # Weight matrix includes weights from tumor types regardless whether edge is significant or not
  mx = match(unionSet,indat[,1]) 
  notNA.ix = which(!is.na(mx))
  notNA = mx[!is.na(mx)]
  
  # Create weight vector and then place in weight matrix
  weightVec = rep(0,length(unionSet))
  weightVec[notNA.ix] = indat[notNA,3] * indat[notNA,4]
  weightMat[,i] = weightVec

  # Binary matrix denotes whether edge passed the CUTOFF 425 in that particular tumor type
  mxB = match(unionSet,indEdgeset[[i]]) 
  signifBinVec = rep(0,length(unionSet))
  signifBinVec[which(!is.na(mxB))] = 1
  signifBinMat[,i] = signifBinVec
  
}#end for i

## SAVE THESE
write.table(weightMat,"PPI_weight_matrix_TOP6_425edges.txt",sep="\t",quote=F)
write.table(signifBinMat,"PPI_binary_matrix_TOP6_425edges.txt",sep="\t",quote=F)
######################################################################


