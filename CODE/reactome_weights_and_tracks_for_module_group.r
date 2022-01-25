########################################################################################
# Average interaction strength for Reactome gene lists
#########################################################################################
# Dependencies for variables
# MAIN.r
#########################################################################################
# REACTOME human pathways
pathways = readRDS("DATA/REACTOME/ReactomePathways_vHuman.rds") # 1669 pathways
# Unique genes
uniq = as.character(read.csv("DATA/REACTOME/unique_genes_Nov11_2015.txt",sep="\t")[,1]) # S8 Table
bp.input = read.csv("DATA/REACTOME/bp_pancan_input.txt",sep="\t",header=F)
multIx = grep("\\|",bp.input[,3])
# Table of edges
master.info = read.csv("DATA/REACTOME/S1A_edge_info.txt",header=T,sep="\t") # S1A Table
master.binary = read.csv("DATA/REACTOME/S1B_binary.txt",header=T,sep="\t") # S1B Table
# Reactome files
reactomeList = read.csv("DATA/REACTOME/ReactomePathwayLabels_vHuman.txt",sep="\t",header=F) 
reactomeHierarchy = read.csv("DATA/REACTOME/ReactomePathwaysRelation_vHuman.txt",sep="\t",header=F)

#########################################################################################
# Genelist by tumor matrix
pathBYtumor = matrix(NA,nrow = length(pathways),ncol=length(types))
colnames(pathBYtumor) = types
rownames(pathBYtumor) = names(pathways)

# Genes corresponding to antibodies
genes = as.character(bp.input[,3])
genes[106:107] = "MTOR" # FRAP1 not gene name
genes[174] = "TIGAR" # C12orf5 not gene name

for(w in 1:length(types)){
  print(types[w])
  # Edges that pass the threshold in tumor type w
  binaries = master.binary[,which(colnames(master.binary)==types[w])]
  ix = which(binaries==1)
  ab.pair = master.binary[ix,1:2]
  weights = master.info[ix,which(colnames(master.info)==types[w])]
  
  # Get corresponding genes for antibody pair
  mx = match(as.character(ab.pair[,1]),prots)
  g1 = as.character(genes[mx])
  mx = match(as.character(ab.pair[,2]),prots)
  g2 = as.character(genes[mx])
  
  bigmat = matrix(NA,ncol=nrow(ab.pair),nrow=length(pathways))
  rownames(bigmat) = names(pathways)
  colnames(bigmat) = paste(ab.pair[,1],ab.pair[,2],sep="&") 
  isFound = bigmat
  isFound[which(is.na(isFound))] = 0
  
  for(i in 1:ncol(bigmat)){
    #print(i)
    genes1 = unlist(strsplit(g1[i],"\\|"))
    genes2 = unlist(strsplit(g2[i],"\\|"))
    
    # Don't allow self interactions
    # Self-interaction is where genes1 and genes2 have a non-empty intersection set
    # i.e there is at least one gene in the interaction set
    if(length(intersect(genes1,genes2))==0){
      for(j in 1:length(pathways)){
        found1 = length(intersect(genes1,pathways[[j]]))
        #if(found1>0){print(j)}
        found2 = length(intersect(genes2,pathways[[j]]))
        #if(found2>0){print(j)}
        if(found1 > 0 & found2 > 0){
          isFound[j,i] = 1
          # penalize by the size of the pathway list
          bigmat[j,i] = abs(weights[i]) / log(length(pathways[[j]])) # absolute weight
        }#end if
      }#end for j
    }# end if
  }#end for i
  
  # Pathway averages (only non-NA entries factor into average)
  p.ave = apply(bigmat,1,mean,na.rm=T)
  pathBYtumor[,w] = p.ave
  # save isFound
  saveRDS(isFound,paste("Reactome_edge_mapping_binary_",types[w],".rds",sep=""))
}#end for w
saveRDS(pathBYtumor,"REACTOME_pathway_by_tumortype_master_matrix_2015Dec10.rds")
###################################################################
#change NANs to 0, and remove 0 rows
X = pathBYtumor
rm(pathBYtumor)
X[which(is.na(X))] = 0

sumX = apply(X,1,sum) # sum
maxX = apply(X,1,max) # max
XX = X[which(maxX > 0),]
dim(XX) # 339 pathways remaining

write.table(XX,"REACTOME_pathway_by_tumortype_N339.txt",sep="\t",quote=F)

###################################################################
# Parent-child analysis
###################################################################

# Up one level
mx = match(child[,2],reactomeHierarchy[,2])
parent.pid =  as.character(reactomeHierarchy[mx,1])
my = match(parent.pid,reactomeList[,1])
parent.path = as.character(reactomeList[my,2])
parent1 = cbind(parent.path,parent.pid)

# Up another level
mx = match(parent1[,2],reactomeHierarchy[,2])
parent.pid =  as.character(reactomeHierarchy[mx,1])
my = match(parent.pid,reactomeList[,1])
parent.path = as.character(reactomeList[my,2])
parent2 = cbind(parent.path,parent.pid)

# Up another level
mx = match(parent2[,2],reactomeHierarchy[,2])
parent.pid =  as.character(reactomeHierarchy[mx,1])
my = match(parent.pid,reactomeList[,1])
parent.path = as.character(reactomeList[my,2])
parent3 = cbind(parent.path,parent.pid)

# Up another level
mx = match(parent3[,2],reactomeHierarchy[,2])
parent.pid =  as.character(reactomeHierarchy[mx,1])
my = match(parent.pid,reactomeList[,1])
parent.path = as.character(reactomeList[my,2])
parent4 = cbind(parent.path,parent.pid)

# Up another level
mx = match(parent4[,2],reactomeHierarchy[,2])
parent.pid =  as.character(reactomeHierarchy[mx,1])
my = match(parent.pid,reactomeList[,1])
parent.path = as.character(reactomeList[my,2])
parent5 = cbind(parent.path,parent.pid)

mx = match(parent5[,2],reactomeHierarchy[,2])
parent.pid =  as.character(reactomeHierarchy[mx,1])
my = match(parent.pid,reactomeList[,1])
parent.path = as.character(reactomeList[my,2])
parent6 = cbind(parent.path,parent.pid)

mx = match(parent6[,2],reactomeHierarchy[,2])
parent.pid =  as.character(reactomeHierarchy[mx,1])
my = match(parent.pid,reactomeList[,1])
parent.path = as.character(reactomeList[my,2])
parent7 = cbind(parent.path,parent.pid)

mx = match(parent7[,2],reactomeHierarchy[,2])
parent.pid =  as.character(reactomeHierarchy[mx,1])
my = match(parent.pid,reactomeList[,1])
parent.path = as.character(reactomeList[my,2])
parent8 = cbind(parent.path,parent.pid)

mx = match(parent8[,2],reactomeHierarchy[,2])
parent.pid =  as.character(reactomeHierarchy[mx,1])
my = match(parent.pid,reactomeList[,1])
parent.path = as.character(reactomeList[my,2])
parent9 = cbind(parent.path,parent.pid)
# PARENT 9 is all NA, so use PARENT 8 as the last parent

# write out table with all parental relationships
B = cbind(CHILD.PID=child[,2],CHILD=child[,1],PARENT1=parent1[,1],PARENT2=parent2[,1],PARENT3=parent3[,1],PARENT4=parent4[,1],
          PARENT5=parent5[,1],PARENT6=parent6[,1],PARENT7=parent7[,1],PARENT8=parent8[,1])
write.table(B,"Edge339_parent_reactome_pathways.txt",sep="\t",quote=F,row.names=F)
######################################################################
# want to add module and group info to heatmap
######################################################################
mxx = match(B[,2],names(pathways))

for(i in 1:length(types)){
  print(types[i])
  isFound = readRDS(paste("Reactome_edge_mapping_binary_",types[i],".rds",sep=""))
  x = isFound[mxx,]
  modmat = grmat = x # result matrix
  
  w = which(master.binary[,which(colnames(master.binary)==types[i])]==1)
  modules = as.character(master.info[w,17])
  groups = as.character(master.info[w,16])
  
  for(j in 1:ncol(x)){
    ww = which(x[,j]==1)
    if(length(ww) > 0){
      modmat[ww,j] = modules[j]
      grmat[ww,j] = groups[j]
    }#end if
  }#end for j
  
  plotmodmat = matrix(0,ncol=6,nrow=nrow(x))
  rownames(plotmodmat) = rownames(x)
  colnames(plotmodmat) = paste("Module",1:6,sep=".")
  
  plotgrmat = matrix(0,ncol=3,nrow=nrow(x))
  rownames(plotgrmat) = rownames(x)
  colnames(plotgrmat) = c("Positive dominant","Negative dominant","Heterogeneous")
  
  for(j in 1:nrow(plotmodmat)){
    tab = table(modmat[j,])
    foundModules = names(tab)
    mx = match(foundModules,1:6)
    plotmodmat[j,mx[!is.na(mx)]] = as.numeric(tab[which(!is.na(mx))])
  }#end for j
  
  for(j in 1:nrow(plotgrmat)){
    tab = table(grmat[j,])
    foundModules = names(tab)
    mx = match(foundModules,colnames(plotgrmat))
    plotgrmat[j,mx[!is.na(mx)]] = tab[which(!is.na(mx))]
  }#end for j
  
  jnk = cbind(GENELIST=rownames(plotmodmat),plotmodmat)
  write.table(jnk,paste("Per_genelist_module_counts_",types[i],".txt",sep=""),sep="\t",quote=F,row.names=F)
  
  jnk = cbind(GENELIST=rownames(plotgrmat),plotgrmat)
  write.table(jnk,paste("Per_genelist_group_counts_",types[i],".txt",sep=""),sep="\t",quote=F,row.names=F)
  
}#end for i


summ = matrix(0,nrow=nrow(B),ncol=ncol(plotmodmat))
for(i in 1:length(types)){
  indat = read.csv(paste("Per_genelist_module_counts_",types[i],".txt",sep=""),sep="\t",header=T)
  summ = summ + indat[,2:ncol(indat)]
}#end for i
mod.ave = summ / length(types)
rownames(mod.ave) = indat[,1]
write.table(mod.ave,"Per_genelist_module_counts_PANCAN_AVERAGE.txt",sep="\t",quote=F)
saveRDS(mod.ave,"../data/Per_genelist_module_counts_PANCAN_AVERAGE.rds")

summ = matrix(0,nrow=nrow(B),ncol=ncol(plotgrmat))
for(i in 1:length(types)){
  indat = read.csv(paste("Per_genelist_group_counts_",types[i],".txt",sep=""),sep="\t",header=T)
  summ = summ + indat[,2:ncol(indat)]
}#end for i
gr.ave = summ / length(types)
rownames(gr.ave) = indat[,1]
write.table(gr.ave,"Per_genelist_group_counts_PANCAN_AVERAGE.txt",sep="\t",quote=F)
saveRDS(gr.ave,"../data/Per_genelist_group_counts_PANCAN_AVERAGE.rds")


####################################################################################################
# Do for all discovery set interactions
####################################################################################################

ab.pair = master.info[,1:2]

# Get corresponding genes for antibody pair
mx = match(as.character(ab.pair[,1]),prots)
g1 = as.character(genes[mx])
mx = match(as.character(ab.pair[,2]),prots)
g2 = as.character(genes[mx])

isFound = matrix(0,ncol=nrow(ab.pair),nrow=length(pathways))
rownames(isFound) = names(pathways)
colnames(isFound) = paste(ab.pair[,1],ab.pair[,2],sep="&") 

for(i in 1:ncol(bigmat)){
  #print(i)
  genes1 = unlist(strsplit(g1[i],"\\|"))
  genes2 = unlist(strsplit(g2[i],"\\|"))
  
  # Don't allow self interactions
  # Self-interaction is where genes1 and genes2 have a non-empty intersection set
  # i.e there is at least one gene in the interaction set
  if(length(intersect(genes1,genes2))==0){
    for(j in 1:length(pathways)){
      found1 = length(intersect(genes1,pathways[[j]]))
      #if(found1>0){print(j)}
      found2 = length(intersect(genes2,pathways[[j]]))
      #if(found2>0){print(j)}
      if(found1 > 0 & found2 > 0){
        isFound[j,i] = 1
      }#end if
    }#end for j
  }# end if
}#end for i

outmat = matrix(NA,nrow=ncol(isFound),ncol=1)
rownames(outmat) = colnames(isFound)

for(i in 1:nrow(outmat)){
  w = which(isFound[,i]==1)
  #print(length(w))
  if(length(w) > 0){
    outmat[i] = paste(rownames(isFound)[w],collapse=", ")
  }#end if
}#end for i

write.table(outmat,"Genelists_corresponding_to_discovery_set_interactions.txt",sep="\t",quote=F,row.names=F)





