##################################################################################################
# PCA, hierarchical clustering, and heatmap for the discovery set 
##################################################################################################
# Dependencies for variables
# MAIN.r, create_discovery_set.r
##################################################################################################
# Principal components on tumor types
# Hierarchical clustering on edges
pc = prcomp(t(weightMat))
hc = hclust(dist(weightMat),method="ward.D2")
classes = cutree(hc,3)
#outfile = cbind(names(classes),classes)
#write.table(outfile,"cluster_info.txt",sep="\t",quote=F,row.names=F)
##################################################################################################
# Heatmap snippet
red = rgb(1,0,0); blue = rgb(0,0,1); white=rgb(1,1,1)
bluered = colorRampPalette(c(blue,white,red))
hclustWard = function(d){
  hclust(d,method="ward.D2")
}#end 
##################################################################################################
# Transpose data
tW = t(weightMat)
recurrence = apply(signifBinMat,1,sum)
##################################################################################################
# Greens for recurrence
green=rgb(0,1,0); whitegreen = colorRampPalette(c(white,green))
greenHues = whitegreen(11)
annotation = matrix(NA,ncol=1,nrow=ncol(tW))
annotation[,1] = greenHues[recurrence]
colnames(annotation)[1] = "RECURRENCE"

tumorOrder = c(2,1,3,4,8,5,7,6,10,9,11)
pdf("Discovery_set_heatmap_ward_linkage.pdf",width=10)
heatmap.3(tW[tumorOrder,],Rowv=FALSE,dendrogram="column",na.color=gray(7/8),
          hclustfun=hclustWard,col=bluered,key=T,keysize=1,margins=c(10,4),
          labCol=NA,cexRow=1.2,ColSideColors=annotation,NumColSideColors = 1,side.height.fraction=0.2)
dev.off()
##################################################################################################








