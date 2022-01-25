########################################################################
# x is the output from communities() in igraph
# prots is the list of all elements in x (only those in x)
convertCommunityToAdjacency = function(x,prots){
  conmat = matrix(0,nrow=length(uniqs),ncol=length(prots))
  rownames(conmat) = colnames(conmat) = prots
  for(i in 1:length(x)){
    this = x[[i]]
    if(length(this) > 1){
      mx = match(this,prots)
      ind1 = t(combn(mx,2))
      ind2 = ind1[,2:1]
      conmat[ind1] = conmat[ind1] + 1 # if communities() works correctly, all communities have to be mutually exclusive so conmat[ind1] has to be 0 all the time (before adding 1)
      conmat[ind2] = conmat[ind2] + 1
    }#end if
   }#end for i
  return(conmat)
}#end for function
########################################################################

getOverlapInTopX = function(unionList,types,topX){
  intset = matrix(NA,length(unionList),length(unionList))
  colnames(intset) = rownames(intset) = types
  for(i in 1:length(types)){
    #print(i)
    thisUnion = unionList[[i]]
    pu = paste(thisUnion[1:topX,1],thisUnion[1:topX,2],sep=".")
    for(j in 1:length(types)){
      bb = unionList[[j]]
      bbp = paste(bb[1:topX,1],bb[1:topX,2],sep=".")
      lenint = length(intersect(pu,bbp))
      lenunion = length(union(pu,bbp))
      intset[i,j] = lenint / lenunion
    }#edn for j
  }#end for
  return(overlap.matrix=intset)
}#end function

########################################################################

longPrecisionRecall = function(edgelist,benchmark,prots,recall.limit=1,sparsePred=FALSE,node1ind=2,node2ind=3,returnEdges=TRUE){
  
  thispred = edgelist[,c(node1ind,node2ind)]
  numProts = length(prots)
  # Is the number of predicted edges the same as the number of all possible edges?
  # If no, then sparsePred is TRUE, so you bind a random component to the end of predictions
  if(sparsePred){
      fullset = t(combn(1:numProts,2))
      fullsetdot =  paste(fullset[,1],fullset[,2],sep=".")
      thispreddot = paste(thispred[,1],thispred[,2],sep=".")
      remaining = setdiff(fullsetdot,thispreddot)
      
      unifs = runif(length(remaining),0,1)
      unifOrder = order(unifs)
      nums = unlist(strsplit(remaining[unifOrder],"\\."))
      node1 = as.numeric(nums[c(TRUE,FALSE)])
      node2 = as.numeric(nums[c(FALSE,TRUE)])
      randomEdges = cbind(node1,node2)
      thispred = rbind(thispred,randomEdges) 
  }#end if
  
  # Ensure alphabetical order in BP_PRIOR output
  ind = which(!(benchmark[,1] < benchmark[,2]))
  benchmark[ind,] = benchmark[ind,c(2,1)]
  # Ensure alphabetical order in computational predictions
  ce = cbind(prots[thispred[,1]],prots[thispred[,2]])
  ind = which(!(ce[,1] < ce[,2]))
  ce[ind,] = ce[ind,c(2,1)]
  #Combine the edges with a dot
  pcdot = paste(benchmark[,1],benchmark[,2],sep=".")
  cedot = paste(ce[,1],ce[,2],sep=".")
  is.match = match(cedot,pcdot)
  is.match[!is.na(is.match)] = 1
  is.match[is.na(is.match)] = 0
  
  cumulative = cumsum(is.match)
  precision = cumulative / (1:length(cumulative))
  recall = cumulative / length(pcdot)
  
  cut.ind = max(which(recall <= recall.limit))
  precision = precision[1:cut.ind]
  recall = recall[1:cut.ind]
  is.match = is.match[1:cut.ind]
  
  #area under curve is sum_from_i=1_to_n-1 (p_i+1 + p_i) * (r_i+1 - r_i) / 2
  aupr = sum((precision[-1] + precision[-length(precision)]) * (recall[-1] - recall[-length(recall)]) / 2)
  
  if(returnEdges){
    return(list(area=aupr,prec=precision,recall=recall,is.match=is.match,sortedEdges=thispred))
  }else{
    return(list(area=aupr,prec=precision,recall=recall,is.match=is.match))
  }#end if-else

}#end function
########################################################################


elasticNetwork = function(X,nfold=10,alpha=0.5,verbose=FALSE){
  
  X = scale(X)
  B = matrix(0, nrow = ncol(X), ncol = ncol(X))
  colnames(B) = 1:ncol(X)
  pcor =  NULL
  n.opt = alpha.opt = rep(NA,ncol(X))
  
  if (verbose == TRUE) {
    cat(paste("Performing local elastic net regressions\n"))
    cat(paste("Vertex no "))
  }
  for(i in 1:ncol(X)){
    if (verbose == TRUE) {
      if ((i/10) == floor(i/10)) {
        cat(paste(i, "..."))
      }
    }
    M = X[,-i]
    p = X[,i]
    fit = elasticNetGridSearch(M,p,nfold=nfold,alpha=alpha)
    coefi = fit$coefficients
    B[i,-i] = coefi
    
    n.opt[i] = fit$nfold.opt
    alpha.opt[i] = fit$alpha.opt
  }#end for i
  pcor = Beta2parcor(B, verbose = verbose)
  cat(paste("\n"))
  
  return(list(pcor=pcor,n.opt=n.opt,alpha.opt=alpha.opt))
  
}#end function


########################################################################



elasticNetGridSearch = function(x,y,nfold=10,alpha=0.5){
  cv.alpha = lambda.opt = matrix(NA,ncol=length(alpha),nrow=length(nfold))
  #colnames(cv.alpha) = alphaVec
  #rownames(cv.alpha) = kvec
  # grid search involving cross-validation splits and alpha values
  x = as.matrix(x)
  y = as.matrix(y)
  for(i in 1:length(nfold)){
    for(j in 1:length(alpha)){
      res <- cv.glmnet(x,y,nfolds=nfold[i],alpha=alpha[j],intercept=FALSE)
      # cross-validation error when lambda is the largest such that it's within 1se of predictive power
      cv.alpha[i,j] = res$cvm[which(res$lambda == res$lambda.1se)]
      lambda.opt[i,j] = res$lambda.1se
    }#end for j
  }#end for i
  opt.ind = which(cv.alpha == min(cv.alpha),arr.ind=TRUE)
  k.opt = nfold[opt.ind[1]]
  a.opt = alpha[opt.ind[2]]
  lambda = lambda.opt[opt.ind[1],opt.ind[2]]
  model = glmnet(x,y,alpha=a.opt,intercept=FALSE)
  pred = predict(model,newx = x, type="coefficients",s=lambda)
  intercept = pred[1]
  coefficients = pred[-1]
  names(coefficients) = rownames(pred)[-1]
  return(list(lambda.opt=lambda,cv.grid=cv.alpha,nfold.opt=k.opt,alpha.opt=a.opt,intercept=intercept,
              coefficients=coefficients))
}#end function

########################################################################

mylars.vLambda = function (X, y,lambda = 0.01, use.Gram = TRUE, normalize = TRUE) 
{
  x <- X
  n <- length(y)
  if (use.Gram == TRUE) {
    type = "covariance"
  }
  if (use.Gram == FALSE) {
    type = "naive"
  }
  globalfit <- glmnet(x, y, family = "gaussian", lambda = lambda, standardize = normalize, type.gaussian = type)
  coefficients = predict(globalfit, type = "coefficients",s = lambda)
  intercept = coefficients[1]
  coefficients = coefficients[-1]
  names(coefficients) = 1:ncol(X)
  object <- list(lambda = lambda,intercept = intercept, coefficients = coefficients)
  invisible(object)
}#end function

########################################################################

lasso.net.vLambda = function (X, lambda = 0.01, use.Gram = FALSE, verbose = FALSE) 
{
  p <- ncol(X)
  X <- scale(X)
  colnames(X) <- 1:p
  B <- matrix(0, nrow = p, ncol = p)
  colnames(B) <- 1:p
  if (verbose == TRUE) {
    cat(paste("Performing local lasso regressions\n"))
    cat(paste("Vertex no "))
  }
  for (i in 1:p) {
    if (verbose == TRUE) {
      if ((i/10) == floor(i/10)) {
        cat(paste(i, "..."))
      }
    }
    noti <- (1:p)[-i]
    yi <- X[, i]
    Xi <- X[, noti]
    dummy <- mylars.vLambda(Xi, yi,lambda = lambda, use.Gram = use.Gram,normalize=TRUE)
    coefi <- dummy$coefficients
    B[i, -i] <- coefi
  }
  pcor <- Beta2parcor(B, verbose = verbose)
  cat(paste("\n"))
  return(list(pcor = pcor,lambda=dummy$lambda))
}#end function

########################################################################
ridge.net.vLambda = function (X, lambda = NULL, countLambda=500, plot.it = FALSE, scale = TRUE, k = 10, verbose = FALSE) 
{
  if (is.null(lambda) == TRUE) {
    ss <- seq(-10, -1, length = countLambda)
    ss <- 10^ss
    n <- nrow(X)
    nn <- n - floor(n/k)
    lambda <- ss * nn * ncol(X)
  }
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X, scale = scale)
  B <- matrix(0, nrow = p, ncol = p)
  lambda.opt <- rep(0, p)
  cat(paste("Performing local ridge regressions\n"))
  cat(paste("Vertex no "))
  for (i in 1:p) {
    if ((i/10) == floor(i/10)) {
      cat(paste(i, "..."))
    }
    noti <- (1:p)[-i]
    yi <- X[, i]
    Xi <- X[, noti]
    r.cv = ridge.cv(Xi, yi, lambda = lambda, scale = scale, 
                    plot.it = plot.it, k = k)
    B[i, -i] = r.cv$coefficients
    lambda.opt[i] = r.cv$lambda.opt
  }
  pcor <- Beta2parcor(B, verbose = verbose)
  return(list(pcor = pcor))
}
########################################################################

friendly.glmnet = function (type=c("ridge","lasso","elasticnet"),X, lambda = NULL, countLambda = 500, scale = TRUE, k = 10, verbose = FALSE, plot.it = FALSE) 
{
  # if lambda is NULL, generate COUNTLAMBDA values of lambda
  if (is.null(lambda) == TRUE) {
    ss <- seq(-10, -1, length = countLambda)
    ss <- 10^ss
    n <- nrow(X)
    nn <- n - floor(n/k)
    lambda <- ss * nn * ncol(X)
  }
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X, scale = scale)
  B <- matrix(0, nrow = p, ncol = p)
  lambda.opt <- rep(0, p)
  cat(paste("Performing local ridge regressions\n"))
  cat(paste("Vertex no "))
  for (i in 1:p) {
    if ((i/10) == floor(i/10)) {
      cat(paste(i, "..."))
    }
    noti <- (1:p)[-i]
    yi <- X[, i]
    Xi <- X[, noti]
    if(type=="ridge"){
      alpha = 0
    }else if(type=="lasso"){
      alpha=1
    }#### else if elastic net, need to optimize alpha
    
    #glmnet wants a decreasing sequence
    lambdarev = rev(lambda)
    cv.r = cv.glmnet(Xi,yi,family="gaussian",lambda=lambdarev,nfolds=k,alpha=alpha,intercept=FALSE)
    B[i, -i] = coef(cv.r)[-1] # remove intercept
    lambda.opt[i] = cv.r$lambda.min
    #if(plot.it){
    #  plot(cv.r)
    #}
  }#end for
  pcor <- Beta2parcor(B, verbose = verbose)
  return(list(pcor = pcor,lambda.opt=lambda.opt))
}#end function


########################################################################
aupr = function(prec,recall){
  
  # Area under Precision-Recall curve
  # many trapezoids
  # area = sigma_from_i_to_(n-1) {p_i + p_(i+1)} * {r_(i+1) - r_i} / 2
  
  na.ind = which(is.na(prec) | is.nan(prec))
  if(length(na.ind) > 0){
    prec = prec[-na.ind]
    recall = recall[-na.ind]
  }#end if
  #if(any(is.nan(prec))) stop("NaN value detected in precision")
  #if(any(is.nan(recall))) stop("NaN value detected in precision")
  
  trapezoid1 = recall[-1] - recall[-length(recall)]
  trapezoid2 = prec[-length(prec)] + prec[-1] 
  area = sum(trapezoid1 * trapezoid2 / 2)
  return(area)
}#end function
########################################################################
extract.network.vY = function (network.all, method = c("prob", "qval", "number","pval"),cutoff = 0.8, verbose = FALSE) 
{
  method = match.arg(method)

  if (method == "prob") {
    if (cutoff < 0 | cutoff > 1) 
      stop("edges: cutoff for prob must be between 0 and 1")
    edges = network.all[network.all$prob > cutoff, ]
  }
  else if (method == "qval") {
    if (cutoff < 0 | cutoff > 1) 
      stop("edges: cutoff for qval must be between 0 and 1")
    edges = network.all[network.all$qval < cutoff, ]
  }
  else if (method == "number") {
    if (cutoff%%1 != 0) 
      stop("edges: cutoff for \"number\" must be integer")
    edges = network.all[1:cutoff, ]
  }
  else if (method == "pval") {
    if (cutoff < 0 | cutoff > 1) 
      stop("edges: cutoff for pval must be between 0 and 1")
    edges = network.all[network.all$pval < cutoff, ]
  }
  
  if (verbose == TRUE) {
    cat("\nSignificant edges: ", nrow(edges), "\n")
    cat("    Corresponding to ", round(nrow(edges)/nrow(network.all), 
                                       4) * 100, "%  of possible edges \n")
  }
  network = edges

  return(network)
}# end function

#########################################################################
get.pcor.XonY = function(X,Y){
  pcor = matrix(NA,nrow=ncol(Y),ncol=ncol(X))
  for(i in seq(ncol(Y))){
    print(i)
    obj = carscore(X,Y[,i],verbose=F)
    pcor[i,] = as.numeric(obj)
  }#end for i
  qvalvec = fdrtool(c(pcor),statistic="correlation")$qval
  qval = matrix(qvalvec,nrow=ncol(Y),ncol=ncol(X))
  rownames(qval) = colnames(Y)
  colnames(qval) = colnames(X)
  return(list(pcor=pcor,qval=qval))
}#end function

get.out.matrix = function(obj,mirnames,protnames,type=c("protONmir","mirONprot"),cutoff=0.05){
  subs = which(obj$qval <= cutoff)
  inds = which(obj$qval <= cutoff, arr.ind=TRUE)
  
  if(type=="mirONprot"){
    out1 = cbind(mirna=mirnames[inds[,2]],
                 prot=protnames[inds[,1]],
                 pcor=round(obj$pcor[subs],5),
                 qval=round(obj$qval[subs],5))
  }else if(type=="protONmir"){
    out1 = cbind(prot=protnames[inds[,2]],
                 mirna=mirnames[inds[,1]],
                 pcor=round(obj$pcor[subs],5),
                 qval=round(obj$qval[subs],5))
  }else{
    stop("type not recognized")
  }
  
  ord = order(as.numeric(out1[,4]))
  out2 = out1[ord,]
  return(out2)
}#end function

compare.XonY.YonX = function(XonY,YonX){
  list1 = paste(XonY[,1],XonY[,2],sep=".")
  list2 = paste(YonX[,2],YonX[,1],sep=".")
  intlist = intersect(list1,list2)
  resmat = matrix(NA,nrow=length(intlist),ncol=6)
  colnames(resmat) = c("mirna","prot","pcor.mirs.on.prot","qval.mirs.on.prot","pcor.prots.on.mir","qval.prots.on.mir")
  resmat[,1:4] = XonY[match(intlist,list1),]
  resmat[,5:6] = YonX[match(intlist,list2),3:4]
  ord1 = order(as.numeric(resmat[,3])*as.numeric(resmat[,5]),decreasing=TRUE) # order according to product of the partial cors
  resmat2 = resmat[ord1,]
  return(resmat2)
}#end function


################################################
# Find number of triangles
################################################
numTriangle = function(edges,prots){
  
  edgelist = cbind(prots[edges[,1]],prots[edges[,2]])
  net = network(edgelist,directed=FALSE,loops=F,multiple=F)
  x = as.matrix.network(net, matrix.type = "adjacency")
  insert = c()
  
  for(i in 1:ncol(x)){
    #print(i)
    dist1n = which(x[,i]==1)
    len = length(dist1n)
    dist2n = vector("list",len)
    for(j in 1:len){
      dist2n[[j]] = which(x[,dist1n[j]]==1)
      
      thisLen = length(dist2n[[j]])
      
      triangleInd = dist2n[[j]][which(x[i,dist2n[[j]]]==1)]
      # Merge i and dist1n[j] and each member of triangleInd
      if(length(triangleInd) > 0){
        for(k in 1:length(triangleInd)){
          sorted = sort(c(i,dist1n[j],triangleInd[k]))
          insert = rbind(insert,paste(sorted,collapse="."))
        }#end for k
      }#end if
      
    }#end for j
    
  }#end for i
  
  uniqueInsert = unique(insert)
  return(nrow(uniqueInsert))
}#end function

################################################
# 
################################################
plot.network.sequence = function(networkList,from,to,prots,colind,colorNums,algoName,versionName,doEdge=TRUE,doNode=TRUE){
  
  for(i in from:to){
    print(i)
    # ALL SAMPLES
    edges = networkList[[i]]
    edgelist = cbind(prots[edges[,1]],prots[edges[,2]])
    edgewid = 30 * abs(edges[,3])
    edgecol = ifelse(sign(edges[,3]) > 0,2,4)
    net = network(edgelist,directed=FALSE,loops=F,multiple=F)
    x = as.matrix.network(net, matrix.type = "adjacency")
    vertexSize = log(apply(x,2,sum)) + 0.5
    vertexColors = colind[colorNums[match(colnames(x),prots)]]
    
    coords = NA
    
    if(doEdge){
      # Edge-coloring
      pdf(paste("edgeColor_noText_",algoName,"_vRPPA_ovca_",versionName,"_cutoff",i,".pdf",sep=""))
      coords = plot.network(net,displaylabels=F,label.pad=0.5,label.cex=0.5,edge.lwd=edgewid,edge.col=edgecol,
                            vertex.col=colors()[230],vertex.border=1,pad=0,vertex.cex=vertexSize)
      dev.off()
      
      pdf(paste("edgeColor_Text_",algoName,"_vRPPA_ovca_",versionName,"_cutoff",i,".pdf",sep=""))
      plot.network(net,coord=coords,displaylabels=T,label.pad=0.5,label.cex=0.7,edge.lwd=edgewid,edge.col=edgecol,
                   vertex.col=colors()[230],vertex.border=0,pad=0,vertex.cex=vertexSize)
      dev.off()
    }#end if
    
    if(doNode){
      # Node-coloring
      pdf(paste("nodeColor_noText_",algoName,"_vRPPA_ovca_",versionName,"_cutoff",i,".pdf",sep=""))
      if(is.na(coords)){
        coords = plot.network(net,coord=coords,displaylabels=F,label.pad=0.5,label.cex=0.5,edge.lwd=edgewid,edge.col="gray",
                              vertex.col=vertexColors,vertex.border=1,pad=0,vertex.cex=vertexSize) 
      }else{
        plot.network(net,coord=coords,displaylabels=F,label.pad=0.5,label.cex=0.5,edge.lwd=edgewid,edge.col="gray",
                     vertex.col=vertexColors,vertex.border=1,pad=0,vertex.cex=vertexSize)
      }
      dev.off()
      
      pdf(paste("nodeColor_Text_",algoName,"_vRPPA_ovca_",versionName,"_cutoff",i,".pdf",sep=""))
      plot.network(net,coord=coords,displaylabels=T,label.pad=0.5,label.cex=0.7,edge.lwd=edgewid,edge.col="gray",
                   vertex.col=vertexColors,vertex.border=0,pad=0,vertex.cex=vertexSize)
      dev.off()
    }#end if
  }#end for
  
}#end function

################################################
# 
################################################
plot.network.custom.color = function(edges,prots,algoName,versionName,i,networkObj=NA,edgeCol=NA,vertexCol=NA,doEdge=FALSE,doNode=FALSE){
  # Edgelist with protein names, and edge widths
  edgelist = cbind(prots[edges[,1]],prots[edges[,2]])
  edgewid = 30 * abs(edges[,3])
  
  # Default network object
  net = networkObj
  if(is.na(net)){
    net = network(edgelist,directed=FALSE,loops=F,multiple=F)
  }
  
  # Adjacency matrix
  x = as.matrix.network(net, matrix.type = "adjacency")
  vertexSize = log(apply(x,2,sum)) + 0.5
  
  # Default edge color
  if(is.na(edgeCol[1])){
    edgeCol = ifelse(sign(edges[,3]) > 0,2,4)
  }
  
  # Default vertex color
  if(is.na(vertexCol[1])){
    vertexCol = colind[colorNums[match(colnames(x),prots)]]
  }
  
  coords = NA
  
  if(doEdge){
    # Edge-coloring
    pdf(paste("edgeColor_noText_",algoName,"_ovca_",versionName,"_cutoff",i,".pdf",sep=""))
    coords = plot.network(net,displaylabels=F,label.pad=0.5,label.cex=0.5,edge.lwd=edgewid,edge.col=edgeCol,
                          vertex.col=colors()[230],vertex.border=1,pad=0,vertex.cex=vertexSize)
    dev.off()
    
    pdf(paste("edgeColor_Text_",algoName,"_ovca_",versionName,"_cutoff",i,".pdf",sep=""))
    plot.network(net,coord=coords,displaylabels=T,label.pad=0.5,label.cex=0.7,edge.lwd=edgewid,edge.col=edgeCol,
                 vertex.col=colors()[230],vertex.border=0,pad=0,vertex.cex=vertexSize)
    dev.off()
  }#end if
  
  if(doNode){
    # Node-coloring
    pdf(paste("nodeColor_noText_",algoName,"_ovca_",versionName,"_cutoff",i,".pdf",sep=""))
    if(is.na(coords)){
      coords = plot.network(net,coord=coords,displaylabels=F,label.pad=0.5,label.cex=0.5,edge.lwd=edgewid,edge.col="gray",
                            vertex.col=vertexCol,vertex.border=1,pad=0,vertex.cex=vertexSize) 
    }else{
      plot.network(net,coord=coords,displaylabels=F,label.pad=0.5,label.cex=0.5,edge.lwd=edgewid,edge.col="gray",
                   vertex.col=vertexCol,vertex.border=1,pad=0,vertex.cex=vertexSize)
    }
    dev.off()
    
    pdf(paste("nodeColor_Text_",algoName,"_ovca_",versionName,"_cutoff",i,".pdf",sep=""))
    plot.network(net,coord=coords,displaylabels=T,label.pad=0.5,label.cex=0.7,edge.lwd=edgewid,edge.col="gray",
                 vertex.col=vertexColors,vertex.border=0,pad=0,vertex.cex=vertexSize)
    dev.off()
  }#end if
}#end function


makeTransparent = function(someColor, alpha=100)
{
  #note: always pass alpha on the 0-255 scale
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

ridge.net.vRPPA = function (X, genelist, lambda = NULL, plot.it = FALSE, scale = TRUE, k = 10, verbose = FALSE) 
{
  if (is.null(lambda) == TRUE) {
    ss <- seq(-10, -1, length = 1000)
    ss <- 10^ss
    n <- nrow(X)
    nn <- n - floor(n/k)
    lambda <- ss * nn * ncol(X)
  }
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X, scale = scale)
  B <- matrix(0, nrow = p, ncol = p)
  lambda.opt <- rep(0, p)
  cat(paste("Performing local ridge regressions\n"))
  cat(paste("Vertex no "))
  for (i in 1:p) {
    if ((i/10) == floor(i/10)) {
      cat(paste(i, "..."))
    }
    removeInd = which(genelist==genelist[i])
    Xi = X[, -removeInd]
    yi = X[, i]
    r.cv = ridge.cv(Xi, yi, lambda = lambda, scale = scale, 
                    plot.it = plot.it, k = k)
    B[i, -removeInd] = r.cv$coefficients
    lambda.opt[i] = r.cv$lambda.opt
  }
  pcor <- Beta2parcor(B, verbose = verbose)
  return(list(pcor = pcor))
}

pls.net.vRPPA = function (X, genelist, scale = TRUE, k = 10, ncomp = 15, verbose = FALSE) 
{
  p = ncol(X)
  n = nrow(X)
  genelist = as.character(genelist)
  if (is.null(ncomp)) {
    ncomp = min(n - 1, ncol(X))
  }
  k = floor(k)
  k = max(1, k)
  if (k > n) {
    cat(paste("k exceeds the number of observations. Leave-one-out is applied.\n"))
    k = n
  }
  B = matrix(0, p, p)
  m = vector(length = p)
  cat(paste("Performing local pls regressions\n"))
  kernel = FALSE
  if (n < (p - 1)) {
    kernel = TRUE
  }
  cat(paste("Vertex no "))
  for (i in 1:p) {
    if ((i/10) == floor(i/10)) {
      cat(paste(i, "..."))
    }
    removeInd = which(genelist==genelist[i])
    Xi = X[, -removeInd]
    yi = X[, i]
    fit = penalized.pls.cv(Xi, yi, scale = scale, k = k,ncomp = ncomp)
    B[i, -removeInd] = fit$coefficients
    m[i] = fit$ncomp.opt
  }
  cat(paste("\n"))
  pcor <- Beta2parcor(B, verbose = verbose)
  return(list(pcor = pcor, m = m))
}# end function 

edge2dist = function(edgelist,prots,method="jac"){
  edges = as.matrix(edgelist[,1:2])
  adjmat = matrix(0,nrow=length(prots),ncol=length(prots))
  adjmat[edges] = 1
  adjmat[edges[,2:1]]=1
  colnames(adjmat) = rownames(adjmat) = prots
  diag(adjmat) = 1
  d = vegdist(t(adjmat),method=method)
  return(d)
}#end function


plotMDS = function(edgelist,prots,filename,colvec=NA){
  
  edges = as.matrix(edgelist[,1:2])
  adjmat = matrix(0,nrow=length(prots),ncol=length(prots))
  adjmat[edges] = 1
  adjmat[edges[,2:1]]=1
  colnames(adjmat) = rownames(adjmat) = prots
  
  #### set the diagonal equal to 1 so that an edge between node X and Y will count for the similarity computation b/w X and Y
  diag(adjmat) = 1
  
  ####
  d = vegdist(t(adjmat),method="jac")
  mds = cmdscale(d,k=3,eig=TRUE)
  
  xmax = max(abs(mds$points[,1]))
  ymax = max(abs(mds$points[,2]))
  
  x.ind = 1:length(prots)
  
  ### Text version
  pdf(filename)
  par(mai=c(0.52,0.52,0.52,0.22))
  par(mfrow=c(2,2))
  
  plot(0,0,type="n",xlab="dimension 1",ylab="dimension 2",xlim=c(-xmax,xmax),ylim=c(-ymax,ymax),col=colind[colorNums],main="x=dim1, y=dim2")
  if(any(is.na(colvec))){
    text(mds$points[x.ind,1],mds$points[x.ind,2],labels=prots,font=2,cex=0.5)
  }else{
    text(mds$points[x.ind,1],mds$points[x.ind,2],labels=prots,font=2,cex=0.5,col=colvec)
  }#end if-else
  
  plot(0,0,type="n",xlab="dimension 1",ylab="dimension 3",xlim=c(-xmax,xmax),ylim=c(-ymax,ymax),col=colind[colorNums],main="x=dim1, y=dim3")
  if(any(is.na(colvec))){
    text(mds$points[x.ind,1],mds$points[x.ind,3],labels=prots,font=2,cex=0.5)
  }else{
    text(mds$points[x.ind,1],mds$points[x.ind,3],labels=prots,font=2,cex=0.5,col=colvec)
  }#end if-else
  
  plot(0,0,type="n",xlab="dimension 2",ylab="dimension 3",xlim=c(-xmax,xmax),ylim=c(-ymax,ymax),col=colind[colorNums],main="x=dim2, y=dim3")
  if(any(is.na(colvec))){
    text(mds$points[x.ind,2],mds$points[x.ind,3],labels=prots,font=2,cex=0.5)
  }else{
    text(mds$points[x.ind,2],mds$points[x.ind,3],labels=prots,font=2,cex=0.5,col=colvec)
  }#end if-else
  dev.off()
}#end for function


plotMDS2 = function(edgelist,prots,filename,mult.ind=NA,labels=FALSE,labels.ind=FALSE,inColor="blue"){
  
  edges = as.matrix(edgelist[,1:2])
  adjmat = matrix(0,nrow=length(prots),ncol=length(prots))
  adjmat[edges] = 1
  adjmat[edges[,2:1]]=1
  colnames(adjmat) = rownames(adjmat) = prots
  
  #### set the diagonal equal to 1 so that an edge between node X and Y will count for the similarity computation b/w X and Y
  diag(adjmat) = 1
  
  ####
  d = vegdist(t(adjmat),method="jac")
  mds = cmdscale(d,k=3,eig=TRUE)
  
  xmax = max(abs(mds$points[,1]))
  ymax = max(abs(mds$points[,2]))
  zmax = max(abs(mds$points[,3]))
  
  x.ind = 1:length(prots)
  
  multTargetCol = rep("gray",length(prots))
  multTargetCol[mult.ind] = inColor
  
  tpos = 1
  tcex = 0.7
  tfont = 1
  
  ### points version
  pdf(filename)
  par(mai=c(0.52,0.52,0.52,0.22))
  par(mfrow=c(2,2))
  
  plot(0,0,type="n",xlab="dimension 1",ylab="dimension 2",xlim=c(-xmax,xmax),ylim=c(-ymax,ymax),col=colind[colorNums],main="x=dim1, y=dim2")
  if(any(is.na(mult.ind))){
    points(mds$points[x.ind,1],mds$points[x.ind,2],pch=16,cex=0.8)
  }else{
    points(mds$points[setdiff(x.ind,mult.ind),1],mds$points[setdiff(x.ind,mult.ind),2],pch=16,cex=0.8,col=multTargetCol[setdiff(x.ind,mult.ind)])
    points(mds$points[mult.ind,1],mds$points[mult.ind,2],pch=1,lwd=1,cex=1,col=inColor)
    if(labels[1] != FALSE){
      text(mds$points[labels.ind,1],mds$points[labels.ind,2],labels=labels,col=inColor,pos=tpos,cex=tcex,font=tfont)
    }
  }#end if-else
  
  plot(0,0,type="n",xlab="dimension 1",ylab="dimension 3",xlim=c(-xmax,xmax),ylim=c(-zmax,zmax),col=colind[colorNums],main="x=dim1, y=dim3")
  if(any(is.na(mult.ind))){
    points(mds$points[x.ind,1],mds$points[x.ind,3],pch=16,cex=0.8)
  }else{
    points(mds$points[setdiff(x.ind,mult.ind),1],mds$points[setdiff(x.ind,mult.ind),3],pch=16,cex=0.8,col=multTargetCol[setdiff(x.ind,mult.ind)])
    points(mds$points[mult.ind,1],mds$points[mult.ind,3],pch=1,lwd=1,cex=1,col=inColor)
    if(labels[1] != FALSE){
      text(mds$points[labels.ind,1],mds$points[labels.ind,3],labels=labels,col=inColor,pos=tpos,cex=tcex,font=tfont)
    }
  }#end if-else
  
  plot(0,0,type="n",xlab="dimension 2",ylab="dimension 3",xlim=c(-ymax,ymax),ylim=c(-zmax,zmax),col=colind[colorNums],main="x=dim2, y=dim3")
  if(any(is.na(mult.ind))){
    points(mds$points[x.ind,2],mds$points[x.ind,3],pch=16,cex=0.8)
  }else{
    points(mds$points[setdiff(x.ind,mult.ind),2],mds$points[setdiff(x.ind,mult.ind),3],pch=16,cex=0.8,col=multTargetCol[setdiff(x.ind,mult.ind)])
    points(mds$points[mult.ind,2],mds$points[mult.ind,3],pch=1,lwd=1,cex=1,col=inColor)
    if(labels[1] != FALSE){
      text(mds$points[labels.ind,2],mds$points[labels.ind,3],labels=labels,col=inColor,pos=tpos,cex=tcex,font=tfont)
    }
  }#end if-else
  dev.off()
}#end for function

getSevenNumPerfCustom = function(pred,trueInt,prots,lenCustomSet){
  num = nrow(pred)
  truth = apply(trueInt,1,paste,collapse=".")
  
  # convert numbers to antibody names
  nameMap = matrix(NA,nrow=num,ncol=2)
  nameMap[,1] = prots[pred[,1]]
  nameMap[,2] = prots[pred[,2]]
  
  predEdges.L2R = apply(nameMap,1,paste,collapse=".")
  # The below line does not work when there is only 1 edge
  if(nrow(nameMap) != 1){
    predEdges.R2L = apply(nameMap[,c(2,1)],1,paste,collapse=".")
  }else{
    predEdges.R2L = paste(rev(nameMap),collapse=".")
  }
  
  # match now
  matched1 = match(truth,predEdges.L2R)
  matched2 = match(truth,predEdges.R2L)
  
  TP = length(which(!is.na(matched1))) + length(which(!is.na(matched2)))
  FP = num - TP
  P = length(truth)
  FN = P - TP
  
  x = lenCustomSet
  allPossible = x * length(prots) - (x*(x+1)/2)
  N = allPossible - P
  TN = N - FP
  
  TPR = TP / P
  FPR = FP / N
  
  PPV = TP / (TP + FP)
  
  return(list(TP=TP,FP=FP,TN=TN,FN=FN,TPR=TPR,FPR=FPR,PPV=PPV))
}#end function

get.random.truth = function(i,tumor,count,prots){
  
  ## Choose random COUNT proteins from PROTS
  rand.ind = sort(sample(1:length(prots),count))
  
  ## PC and BioGRID benchmark
  load(paste("../Feb4_2014_prec_recall_merge_best/truth/ints.",tumor[i],".PC2mergeWithGenemaniaPhysical.rdata",sep="")) # merged
  truth = cbind(match(merged[,1],prots),match(merged[,2],prots))
  # Get only the desired group of antibodies from TRUTH
  newTruthInd = truth[sort(union(which(!is.na(match(truth[,1],rand.ind))),which(!is.na(match(truth[,2],rand.ind))))),]
  newTruth = cbind(prots[newTruthInd[,1]],prots[newTruthInd[,2]])
  
  return(list(edgelist=newTruth,randInd=rand.ind))
  
}#End function get.random.truth


get.RTK.truth = function(i,tumor,prots){
  
  ## RTKs: EGFR, HER2, HER3s, IGF-1R, c-Met, VEGFR2, c-Kit
  rtk.ind = c(grep("EGFR-",prots),grep("EGFR_",prots),grep("HER2-",prots),grep("HER2_",prots),
              grep("HER3-",prots),grep("HER3_",prots),grep("IGF-1R-",prots),grep("IGF-1R_",prots),
              grep("c-Met-",prots),grep("c-Met_",prots),grep("c-Kit-",prots),grep("c-Kit_",prots),
              grep("VEGFR2-",prots),grep("VEGFR2_",prots))
  rtks = matrix(NA,ncol=2,nrow=length(rtk.ind))
  rtks[,1] = rtk.ind
  rtks[,2] = ifelse(grepl("pY",prots[rtk.ind]),"p","t")
  colnames(rtks) = c("ind","type")
  rownames(rtks) = prots[rtk.ind]
  
  ## PC and BioGRID benchmark
  load(paste("../Feb4_2014_prec_recall_merge_best/truth/ints.",tumor[i],".PC2mergeWithGenemaniaPhysical.rdata",sep="")) # merged
  truth = cbind(match(merged[,1],prots),match(merged[,2],prots))
  # Get only the desired group of antibodies from TRUTH
  newTruthInd = truth[sort(union(which(!is.na(match(truth[,1],rtks[,1]))),which(!is.na(match(truth[,2],rtks[,1]))))),]
  newTruth = cbind(prots[newTruthInd[,1]],prots[newTruthInd[,2]])
  
  return(list(edgelist=newTruth,rtks=rtks))
  
}#End function get.RTK.truth



getEdgeList = function(dat,ind,K=0.8,filename=NULL){
  datsub= dat[ind,]
  ggm = ggm.estimate.pcor(datsub)
  perf =  performance.pcor.vY(ggm,cutoff.ggm=K)
  if(!is.null(filename)){
    a = writeConnectivityMatrix(perf,filename)
    return(list(outfile=a,perf=perf))
  }else{
    return(list(outfile=NULL,perf=perf))
  }
}#end function


writeConnectivityMatrix = function(perfObj,filename){
  net = perfObj$net
  connect = perfObj$connectivity
  thisConnect = matrix(connect[unlist(c(net[,2:3]))],ncol=2)
  colnames(thisConnect) = c("degree1","degree2")
  outfile = cbind(net,thisConnect)
  outfile[,2:3] = matrix(rownames(ggm.AB)[unlist(c(net[,2:3]))],ncol=2)
  outfile2 = outfile[,c(2,3,7,8,1,4,5,6)]
  write.table(outfile2,file=filename,sep="\t",quote=F,row.names=F)
  return(outfile=outfile2)
}#end function

writeConnectivityMatrix2 = function(ggm,perfObj,filename){
  net = perfObj$net
  connect = perfObj$connectivity
  thisConnect = matrix(connect[unlist(c(net[,2:3]))],ncol=2)
  colnames(thisConnect) = c("degree1","degree2")
  outfile = cbind(net,thisConnect)
  outfile[,2:3] = matrix(rownames(ggm)[unlist(c(net[,2:3]))],ncol=2)
  outfile2 = outfile[,c(2,3,7,8,1,4,5,6)]
  write.table(outfile2,file=filename,sep="\t",quote=F,row.names=F)
  return(outfile=outfile2)
}#end function

maf2matrix.vRppa = function(MUT,prot.samps,verbose=FALSE){
  
  # Make mutation matrix using RPPA samples
  
  uniq.mut.genes = unique(MUT[,1])
  num.mut.gene = length(uniq.mut.genes)
  mat = matrix(NA,nrow=length(prot.samps),ncol=num.mut.gene)
  
  rownames(mat) = prot.samps
  colnames(mat) = uniq.mut.genes
  
  for(m in 1:num.mut.gene){
    genename = as.character(uniq.mut.genes[m])
    if(verbose==TRUE){
      print(genename)
    }#end if  
    # Can have more than one mutation per sample, so get unique sample names
    mut.samp = unique(MUT[which(MUT[,1]==genename),7])
    # mutStat a binary vector showing mutation status where the
    # samples are ordered according to RPPA samples
    mutStat = match(prot.samps,mut.samp)
    mutStat[which(!is.na(mutStat))] = 1
    mutStat[which(is.na(mutStat))] = 0
    
    mat[,m] = mutStat
    
  }#end for
  
  return(mat)
  
}#end function

formatFirehoseRnaseq = function(dat){
  mat = matrix(as.numeric(as.matrix(dat[-1,-1])),nrow=nrow(dat)-1,ncol=ncol(dat)-1)
  # Get gene names
  vec = unlist(strsplit(as.character(dat[,1]),"\\|"))
  genes = vec[c(FALSE,TRUE)]
  rownames(mat) = genes
  # Get colnames
  colnames(mat) = colnames(dat)[-1]
  return(mat)
}#end function


format.maf = function(mutdat,sortBy=c("sample","gene"),datatype="single"){
  accepted.mut = c("Missense_Mutation","Nonsense_Mutation","Splice_Site","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Translation_Start_Site")
  mx = match(mutdat[,9],accepted.mut)
  remove1 = which(is.na(mx))
  remove2 = which(mutdat[,1]=="Unknown")
  if(datatype=="single"){
    mutnew = mutdat[-union(remove1,remove2),c(1,2,5,6,7,8,9,16,42)]
  }else if(datatype=="pancan"){
    mutnew = mutdat[-union(remove1,remove2),c(1,2,5,6,7,8,9,16,48)]
  }else if(datatype=="DCC"){
    mutnew = mutdat[-union(remove1,remove2),c(1,2,5,6,7,8,9,16)]
  }else{
    stop("Can't recognize datatype variable")
  }#end if-else
  mutnew[,8] = gsub("-",".",substr(mutnew[,8],1,12))
  
  if(sortBy=="sample"){
    mutnew2 = mutnew[order(mutnew[,8]),]
    rownames(mutnew2) = 1:nrow(mutnew2) 
  }else if(sortBy=="gene"){
    mutnew2 = mutnew[order(mutnew[,1]),]
    rownames(mutnew2) = 1:nrow(mutnew2)  
  }else{
    stop("Can't recognize sortBy variable")
  }#end if-else
  return(mutnew2) 
  
}#end function



format.rppa = function(indat){

  # Takes Firehose-like rppa data as input, and converts it to a cbio-portal-like format
  
  # remove whitespace from gene names
  names = gsub(" \\|","|",indat[,1])
  a = strsplit(as.character(names),"\\|")
  type = gene = antibody = rep(NA,length(a))
  for(i in 1:length(a)){
    gene[i] = a[[i]][1]
    antibody[i] = a[[i]][2]
    pS = length(grep("pS",antibody[i]))
    pT = length(grep("pT",antibody[i]))
    pY = length(grep("pY",antibody[i]))
    if((pS > 0 || pT>0) || pY>0){
      type[i] = "phosphoprotein"
    }else{
      type[i] = "protein"
    }#end if-else
  }#end for
  
  # Add slash in between gene names for the antibodies targeting multiple genes
  gene2 = gsub(" ","/",gene)
  mat = cbind(gene=gene2,antibody=antibody,type=type,indat[,2:ncol(indat)])
  
  return(mat)

}#end function


rppaFormat = function(mat){
  mat[which(mat=="GSK3A")] = "GSK3A|GSK3B"
  mat[which(mat=="GSK3B")] = "GSK3A|GSK3B"
  mat[which(mat=="AKT1")] = "AKT1|AKT2|AKT3"
  mat[which(mat=="AKT2")] = "AKT1|AKT2|AKT3"
  mat[which(mat=="AKT3")] = "AKT1|AKT2|AKT3"
  mat[which(mat=="MTOR")] = "FRAP1|MTOR"
  mat[which(mat=="MS4A1")] = "CD20|MS4A1"
  mat[which(mat=="SERPINE1")] = "PAI-1"
  mat[which(mat=="C12orf5")] = "C12ORF5"
  return(mat)
}

getGenemaniaEdges = function(hugo,genemania,hugo2ens){  
  # Given a list of genes, get the genemania interactions from the database variable
  
  # convert hugo to ENSEMBL symbols
  ens = as.character(hugo2ens[match(hugo,hugo2ens[,1]),3])
  ens = ens[-which(ens=="")]
  ens = ens[-which(is.na(ens))]
  
  # 1) Get from GENEMANIA only those rows with entries from ens
  col1match = which(!is.na(match(genemania[,1],ens)))
  col2match = which(!is.na(match(genemania[,2],ens)))
  edge.match = genemania[intersect(col1match,col2match),]
  
  # 2) REmove duplicates from edge.match, also make edges undirected
  edge.match2 = rbind(as.matrix(edge.match),as.matrix(edge.match[,c(2,1,3)]))
  dotted = apply(edge.match2[,1:2],1,paste,collapse=".")
  s.u.d = sort(unique(dotted))
  
  sud2 = c()
  for(i in 1:length(s.u.d)){
    num1 = as.numeric(substr(s.u.d[i],10,15))
    num2 = as.numeric(substr(s.u.d[i],26,31))
    # cat(num1,"and",num2,"\n")
    #print(num2)
    if(num1 < num2){
      sud2 = c(sud2,s.u.d[i]) 
    }#end if
  }#end for
  
  sud2s = strsplit(sud2,"\\.")
  back2hugo = matrix(NA,nrow=length(sud2s),ncol=2)
  for(i in 1:length(sud2)){
    back2hugo[i,1] = as.character(hugo2ens[which(hugo2ens[,3]==sud2s[[i]][1]),1])
    back2hugo[i,2] = as.character(hugo2ens[which(hugo2ens[,3]==sud2s[[i]][2]),1])
  }#end for i
  
  return(back2hugo)
}#end function


getTrueFalseVector = function(pred,trueInt,prots){
  num = nrow(pred)
  truth = apply(trueInt,1,paste,collapse=".")

  # convert numbers to antibody names
  nameMap = matrix(NA,nrow=num,ncol=2)
  nameMap[,1] = prots[pred[,1]]
  nameMap[,2] = prots[pred[,2]]
  
  predEdges.L2R = apply(nameMap,1,paste,collapse=".")
  predEdges.R2L = apply(nameMap[,c(2,1)],1,paste,collapse=".")
  
  # match now
  matched1 = match(predEdges.L2R,truth)
  matched2 = match(predEdges.R2L,truth)

  tfvec = rep(FALSE,length(matched1))
  tfvec[which(!is.na(matched1))] = TRUE
  tfvec[which(!is.na(matched2))] = TRUE

  return(tfvec)
}#end function

getSevenNumPerf = function(pred,trueInt,prots){
  num = nrow(pred)
  #print(num)
  truth = apply(trueInt,1,paste,collapse=".")
  
  # convert numbers to antibody names
  nameMap = matrix(NA,nrow=num,ncol=2)
  nameMap[,1] = prots[pred[,1]]
  nameMap[,2] = prots[pred[,2]]
  
  predEdges.L2R = apply(nameMap,1,paste,collapse=".")
  # The below line does not work when there is only 1 edge
  if(nrow(nameMap) != 1){
    predEdges.R2L = apply(nameMap[,c(2,1)],1,paste,collapse=".")
  }else{
    predEdges.R2L = paste(rev(nameMap),collapse=".")
  }
  
  # match now
  matched1 = match(truth,predEdges.L2R)
  matched2 = match(truth,predEdges.R2L)
  
  TP = length(which(!is.na(matched1))) + length(which(!is.na(matched2)))
  FP = num - TP
  P = length(truth)
  FN = P - TP
  
  numAB = length(prots)
  N = (numAB * (numAB-1) / 2) - P
  TN = N - FP
  
  TPR = TP / P
  FPR = FP / N
  
  PPV = TP / (TP + FP)
  
  return(list(TP=TP,FP=FP,TN=TN,FN=FN,TPR=TPR,FPR=FPR,PPV=PPV))
}#end function

################################################################################################
performance.pcor.vY = function (inferred.pcor, true.pcor = NULL, fdr = TRUE, cutoff.ggm = 0.8, vSparse=FALSE,
                                verbose = FALSE, plot.it = FALSE) 
{
  p <- ncol(inferred.pcor)
  if (fdr) {
      test.results <- ggm.test.edges(inferred.pcor, verbose = verbose, plot = plot.it)
      prob.results <- diag(p)
      for (i in 1:length(test.results$prob)) {
        prob.results[test.results$node1[i], test.results$node2[i]] <- test.results$prob[i]
        prob.results[test.results$node2[i], test.results$node1[i]] <- test.results$prob[i]
        #dim(prob.results)
      }#end for i
      
      # vSparse=TRUE means "don't look at FDR, get the top X percent of the edges, where X is (1-cutoff.ggm)*100
      if(vSparse){
        numEdge = length(which(abs(test.results[,1]) > 0))
        ## Here, in the case of vSparse, cutoff.ggm means that you want the top (1-cutoff.ggm)*100% of the edges 
        cutoff.number = round((1-cutoff.ggm) * numEdge)
        if(cutoff.number > 0){
          net <- extract.network(test.results, method.ggm="number",cutoff.ggm = cutoff.number)
        }else{
          net <- test.results[NULL,]
        }
      }else{
        net <- extract.network(test.results,cutoff.ggm = cutoff.ggm) 
      }#end if-else
      
      num.selected = nrow(net)
      
      # adjacency matrix
      adj <- diag(p)
      if (num.selected != 0) {
        for (j in 1:num.selected) {
          adj[net[j, 2], net[j, 3]] <- 1
          adj[net[j, 3], net[j, 2]] <- 1
        }#end for
      }#end if 
      
  }else{
    
      # adjacency matrix
      adj <- inferred.pcor
      adj[adj != 0] <- 1
      
      # Edge-list
      symadj = inferred.pcor
      symadj[lower.tri(symadj,diag=TRUE)]=0
      res = cbind(pcor=symadj[symadj != 0],which(symadj != 0,arr.ind=TRUE))
      net = res[order(abs(res[,3]),decreasing=TRUE),]
      
      num.selected = nrow(net)
  }#end if-else fdr

  # How many neighbors does each node have, i.e connectivity
  connectivity <- apply(adj, 2, FUN = function(x) return(sum(x != 0))) - 1
  
  #num.selected <- (sum(abs(adj) > 0) - p)/2
  
  ## Percent of positive edges in the set of all selected edges
  positive.cor <- NULL
  if (num.selected > 0) {
      positive.cor <- sum(net[,1] > 0)/num.selected
  }
  
  num.true <- power <- ppv <- NULL
  if (!is.null(true.pcor)) {
    num.true <- (sum(abs(true.pcor) > 0) - p)/2
    num.true.positives <- (sum(abs(true.pcor * adj) > 0) - 
                             p)/2
    num.false.positives <- num.selected - num.true.positives
    power <- -Inf
    if (num.true > 0) {
      power <- num.true.positives/num.true
    }
    ppv <- -Inf
    if (num.selected > 0) {
      ppv <- num.true.positives/num.selected
    }
  }
  
  auc <- NA
  tpr <- fpr <- NA
  TPR <- FPR <- NA
  if (fdr & (!is.null(true.pcor))) {
    xx <- seq(0, 1, length = 500)
    true.binary = (abs(true.pcor) > 0)
    predicted <- sym2vec(prob.results)
    if (var(predicted) > 0) {
      if (plot.it == TRUE) {
        plot.roc = "ROC"
      }
      if (plot.it == FALSE) {
        plot.roc = NA
      }
      myroc <- ROC(predicted, sym2vec(true.binary), plot = plot.roc)
      auc <- myroc$AUC
      TPR <- myroc$res[, 1]
      FPR <- 1 - myroc$res[, 2]
    }
  }
  
  if (!is.null(true.pcor)) {
    tpr <- power
    if (num.true > 0) {
      fpr <- num.false.positives/((p^2 - p)/2 - num.true)
    }
  }
  
  return(list(num.selected = num.selected, power = power, TPR = TPR, 
              FPR = FPR, tpr = tpr, fpr = fpr, ppv = ppv, adj = adj, 
              connectivity = connectivity, positive.cor = positive.cor, 
              auc = auc,net=net))
}#end function
####################################################################################
#### Convert inverse covariance matrix to partial correlation matrix
invcov2parcor = function(mat){
  vars = diag(mat)
  len = length(vars)
  tmp = pcor = matrix(NA,nrow=len,ncol=len)
  for(m in 1:len){
    tmp[m,] = mat[m,]/sqrt(vars[m])   
  }#end for m
      for(m in 1:len){
        pcor[,m] = tmp[,m]/sqrt(vars[m])   
      }#end for m
          return(pcor)
}#end function
####################################################################################

getTrueFalseStats=function(edges,genes,database.edges,hugo2ens,sep=NULL,allowStateInteractions=FALSE){
  
  TP = FP = rep(NA,nrow(edges))
  
  for(i in 1:nrow(edges)){
    
    cat("i =",i,"\n")
    A = genes[edges[i,1]]
    B = genes[edges[i,2]]
    
    # are there multiple gene names in one entry
    if(!is.null(sep)){
      A = strsplit(A,sep)[[1]]
      B = strsplit(B,sep)[[1]]
    }#end if
    
    # Populate distance array, 0 means interaction of two phosphorylation states
    distarray = c()
    for(j in 1:length(A)){
      for(k in 1:length(B)){
        res = edgeDist(A[j],B[k],database.edges=database.edges,database.type="genemania",hugo2ens=hugo2ens,allowStateInteractions=allowStateInteractions)
        distarray[(j-1)*length(B)+k] = res$d
      }#
    }#end for
    
    if(is.element(FALSE,is.na(distarray))){
      TP[i] = 1
      FP[i] = 0
    }else{
      TP[i] = 0
      FP[i] = 1
    }#end if-else
    
  }#end for
  
  return(list(TP=TP,FP=FP,TN=TN,FN=FN))
  
}#end function


edgeDist = function(geneA,geneB,database.edges=NULL,database.type=c("genemania","pc"),hugo2ens=NULL,allowStateInteractions=TRUE){
  # if database.type is genemania, gene names will be converted to ENSEMBL IDs.
  if(geneA==geneB){
    if(allowStateInteractions){
      return(list(d=0,w=NA)) 
    }else{
      return(list(d=NA,w=NA)) 
    }#end
  }#end 
  
  if(database.type == "genemania"){
    geneA = as.character(hugo2ens[which(hugo2ens[,1]== geneA),3])
    geneB = as.character(hugo2ens[which(hugo2ens[,1]== geneB),3])
  }#end if
  
  tmp = getShortestDist(geneA,geneB,database.edges,hugo2ens,type="single",getMid=FALSE)
  return(list(d=tmp$sdist,w=tmp$weight))
}#end function

getEdges = function(symadj){
  # Function devised for mutual information algorithms like aracne, clr, mrnet
  # symadj is symmetric adjacency matrix where the cells show the strength of the relationship  
  symadj[lower.tri(symadj,diag=TRUE)]=0
  res = cbind(which(symadj >0,arr.ind=TRUE),value=symadj[symadj > 0])
  int = res[order(abs(res[,3]),decreasing=TRUE),]
  return(int)
}#end function

getShortestDistSet = function(a,b,edges,hugo2ens){
	intPairs = cbind(rep(a,each=length(b)),rep(b,length(a)))
	sdist = weight = rep(NA,nrow(intPairs))
	for(j in 1:nrow(intPairs)){
		N = getShortestDist(intPairs[j,1],intPairs[j,2],edges,hugo2ens,type="single")
		sdist[j] = N$sdist
		weight[j] = N$weight
		mid[j] = N$mid
	}#end for
	
	if(any(!is.na(sdist))){
		ind = which.min(sdist)
		return(list(sdist=sdist[ind],weight=weight[ind],mid=mid[ind]))
	}else{
		return(list(sdist=NA,weight=NA,mid=NA))
	}#end if-else
	
}#end function

getShortestDist = function(a,b,edges,hugo2ens=NULL,type="single",getMid=FALSE){
	
	if(type=="set"){
		P = getShortestDistSet(a,b,edges,hugo2ens)
		sdist = P$sdist
		weight = P$weight
		middleNode = P$mid
	}else if(type=="single"){
    
	  if(a==b){
        if(getMid){
    	    return(list(sdist=0,weight=NA,mid=NA))
    	  }else{
    	    return(list(sdist=0,weight=NA))
    	  }#end if
	  }#end if
		
		# Check to see if a and b are found in the same row	
		a.in.col1 = which(edges[,1]==a)	
		b.in.col2 = which(edges[,2]==b)
		row = intersect(a.in.col1,b.in.col2)
		if(length(row) == 0){
			a.in.col2 = which(edges[,2]==a)	
			b.in.col1 = which(edges[,1]==b)
			row = intersect(a.in.col2,b.in.col1)
		}#end if
		
		# If found in the same row, assign sdist as 1, else check if sdist is 2
		if(length(row) == 1){
			sdist = 1
			weight = input[row,3]
			middleNode = NA
		}else{
			#checkIfDist2
			a.neighbors = union(edges[a.in.col1,2],edges[a.in.col2,1])
			b.neighbors = union(edges[b.in.col1,2],edges[b.in.col2,1])
			
			int.neighbors = intersect(a.neighbors,b.neighbors)
			
			if(length(int.neighbors) > 0){
				sdist = 2
				weight = NA
				tmp = match(int.neighbors,hugo2ens[,3])
        if(getMid){
				    middleNode = paste(sort(hugo2ens[tmp[!is.na(tmp)],1]),collapse=", ")
        }#end if
			}else{
				sdist = NA
				weight = NA
				middleNode = NA
			}#end if-else
		}#end if-else
		
	}else{
		stop("shortest distance type not valid")
	}#end if-else

	if(getMid){
    return(list(sdist=sdist,weight=weight,mid=middleNode))
	}else{
	  return(list(sdist=sdist,weight=weight))
	}
}#end function


getEnsembleIDs = function(edges,hugo2ens){

	# some cleaning necessary
	annot1 = as.character(edges[,2])
	annot2 = as.character(edges[,4])
	
	# Remove p from the beginning of phosphoproteins
	ind1 = which(substr(annot1,1,1)=="p")
	annot1[ind1] = substring(annot1[ind1],2)
	ind2 = which(substr(annot2,1,1)=="p")
	annot2[ind2] = substring(annot2[ind2],2)
	
	# Change FRAP1 to MTOR
	annot1[which(annot1 == "FRAP1")] = "MTOR"
	annot2[which(annot2 == "FRAP1")] = "MTOR"
	
	# Change CD20 to MS4A1
	annot1[which(annot1 == "CD20")] = "MS4A1"
	annot2[which(annot2 == "CD20")] = "MS4A1"
	
	# Initialize ENSEMBLE sets
	set1 = vector("list",length(annot1))
	set2 = vector("list",length(annot2))
	
	# Fill in the ones that have only one gene name in the ANNOT entry
	ind3 = grep("/",annot1)
	setdiff1 = setdiff(1:length(annot1),ind3)
	set1[setdiff1] = as.character(hugo2ens[match(annot1[setdiff1],hugo2ens[,1]),3])
	
	ind4 = grep("/",annot2)
	setdiff2 = setdiff(1:length(annot2),ind4)
	set2[setdiff2] = as.character(hugo2ens[match(annot2[setdiff2],hugo2ens[,1]),3])
	
	# Fill in the ones with multiple entries
	vec1 = strsplit(annot1[ind3],"/")
	for(i in 1:length(ind3)){
		set1[[ind3[i]]] = as.character(hugo2ens[match(vec1[[i]],hugo2ens[,1]),3])
	}#end for i
	
	vec2 = strsplit(annot2[ind4],"/")
	for(i in 1:length(ind4)){
		set2[[ind4[i]]] = as.character(hugo2ens[match(vec2[[i]],hugo2ens[,1]),3])
	}#end for i
	
	return(list(set1=set1,set2=set2,gene1=annot1,gene2=annot2))
}#end function



getRppaEdaPlots = function(string,mat1,mat2,redblue){
  
  joint = rbind(mat1,mat2) 
  rownames(joint) = c(rownames(mat1),rownames(mat2))
  # Exploratory Data Analysis
  pdf(paste(string,"_eda1.pdf",sep=""))
  plot(density(as.numeric(mat1[1,])),ylim=c(0,1),xlim=c(-8,8),main="Individual antibody densities",lty=3)
  for(i in 1:nrow(mat1)){
    lines(density(as.numeric(mat1[i,])),lty=3) 
  }#end for
  for(i in 1:nrow(mat2)){
    lines(density(as.numeric(mat2[i,])),col=2,lty=3) 
  }#end for
  legend("topright",legend=c("protein_level","phosphoprotein"),col=c(1,2),lty=1)
  dev.off()
  
  # Heatmaps
  pdf(paste(string,"_protlevel_heatmap.pdf",sep=""))
  heatmap(mat1,labCol=F,col=redblue,xlab="Samples",ylab="Genes")
  dev.off()
  pdf(paste(string,"_phospho_heatmap.pdf",sep=""))
  heatmap(mat2,labCol=F,col=redblue,xlab="Samples",ylab="Genes")
  dev.off()
  pdf(paste(string,"_allantibodies_heatmap.pdf",sep=""))
  heatmap(rbind(mat1,mat2),labCol=F,col=redblue,xlab="Samples",ylab="Genes")
  dev.off()
  
  # PCA
  pc1 = prcomp(mat1)
  plotPCA123(pc1$x,paste(string,"_protlevel_pca",sep=""))
  
  pc2 = prcomp(mat2)
  plotPCA123(pc2$x,paste(string,"_phospho_pca",sep=""))
  
  pc3 = prcomp(joint)
  colors = c(rep(2,nrow(mat1)),rep(3,nrow(mat2)))
  plotPCA123(pc3$x,paste(string,"_allantibodies_pca",sep=""),color=colors)
  
}#end function

plotDecoyPartialCor = function(mat){
	# ncol(mat) - nrow(mat) vectors introduce computational singularity
	# So we add  ncol(mat) - nrow(mat) -1 random vectors as decoy

	numRows2add = ncol(mat) - nrow(mat) -1
	decoy = matrix(rnorm(numRows2add*ncol(mat), mean=0, sd=1),ncol=ncol(mat))
	mat2 = rbind(as.matrix(mat),decoy)
	covMat2 = cor(t(mat2))
	inv2 = solve(covMat2)

	# divide inv2 into 3 sections: partial correlations among antibodies, between antibodies
	# and decoys, and among decoys

	parCor = vector("list",3)
	parCor[[1]] = inv2[1:nrow(mat),1:nrow(mat)] # ab only
	parCor[[2]] = inv2[(nrow(mat)+1):ncol(inv2),1:nrow(mat)] # ab and decoy
	parCor[[3]] = inv2[(nrow(mat)+1):ncol(inv2),(nrow(mat)+1):ncol(inv2)] # decoy only

	abOnly = parCor[[1]][lower.tri(parCor[[1]],diag=T)]
	abDecoy = parCor[[2]][lower.tri(parCor[[2]],diag=T)]
	DecoyOnly = parCor[[3]][lower.tri(parCor[[3]],diag=T)]

	s.DecoyOnly = sort(DecoyOnly)
	s.abDecoy = sort(abDecoy)
	s.abOnly = sort(abOnly)

	xxlim = c(min(inv2),max(inv2))
	yylim = c(0,max(length(abOnly),length(abDecoy),length(DecoyOnly)))
	plot(s.DecoyOnly,1:length(DecoyOnly),col=colors()[536],pch=0,xlim=xxlim,ylim=yylim,xlab="partial correlations",ylab="index",main=paste("[",tumorType[i]," + decoy] partial correlations",sep=""))
	points(s.abDecoy,1:length(abDecoy),col=colors()[229],pch=2)
	points(s.abOnly,1:length(abOnly),col=1)

	maxRandom = max(max(DecoyOnly),max(abDecoy))
	minRandom = min(min(DecoyOnly),min(abDecoy))

	points(s.abOnly[s.abOnly > maxRandom],which(s.abOnly > maxRandom),col=3)
	points(s.abOnly[s.abOnly < minRandom],which(s.abOnly < minRandom),col=3)

	abline(v=maxRandom,col=2,lwd=2)
	abline(v=minRandom,col=2,lwd=2)

	legend("topright",legend=c("within tumor sig.","within tumor non-sig.","tumor vs decoy","decoy vs decoy"),col=c(3,1,colors()[229],colors()[536]),pch=c(1,1,2,0))
}#end function



# plotPCA123.r plots PC1~PC2, PC1~PC3, PC2~PC3 on the same plot
# x is the x object from PCA
# pdfname should not have any extensions at the end
# color is a vector showing color assignments for all samples
# coreVec is a binary vector -same length as color- that has 1 for the core samples
#    and 0 for noncore samples.

plotPCA123 = function(x,pdfname,color=rep(1,numSamp),coreVec=rep(1,numSamp)){
    
    numSamp = nrow(x)   
    
    varSum = sum(apply(x,2,var))
    var1 = round(100 * var(x[,1]) / varSum,digits=2)
    var2 = round(100 * var(x[,2]) / varSum,digits=2)
    var3 = round(100 * var(x[,3]) / varSum,digits=2)

    # Plot the noncore ones first
    noncore = which(coreVec == 0)
    core = which(coreVec == 1)
    myXLIM = c(min(x[,1]),max(x[,1]))
    myYLIM = c(min(x[,2]),max(x[,2]))

    pdf(paste(pdfname,".pdf",sep=""))
    par(mfrow = c(2,2))

    plot(x[noncore,1],x[noncore,2],xlab = paste("PC1 (",var1,"%)",sep=""),ylab=paste("PC2 (",var2,"%)",sep=""),main=paste("PC1 vs PC2, ",pdfname,sep=""),col = color[noncore],xlim=myXLIM,ylim=myYLIM)
    points(x[core,1],x[core,2],pch=19,col=color[core])

    myXLIM = c(min(x[,1]),max(x[,1]))
    myYLIM = c(min(x[,3]),max(x[,3]))
    plot(x[noncore,1],x[noncore,3],xlab = paste("PC1 (",var1,"%)",sep=""),ylab=paste("PC3 (",var3,"%)",sep=""),main=paste("PC1 vs PC3, ",pdfname,sep=""),col = color[noncore],xlim=myXLIM,ylim=myYLIM)
    points(x[core,1],x[core,3],pch=19,col=color[core])

    myXLIM = c(min(x[,2]),max(x[,2]))
    myYLIM = c(min(x[,3]),max(x[,3]))
    plot(x[noncore,2],x[noncore,3],xlab = paste("PC2 (",var2,"%)",sep=""),ylab=paste("PC3 (",var3,"%)",sep=""),main=paste("PC2 vs PC3, ",pdfname,sep=""),col = color[noncore],xlim=myXLIM,ylim=myYLIM)
    points(x[core,2],x[core,3],pch=19,col=color[core])

    dev.off()

}#end function



getShortestDistance = function(pci,int2,ab){
  # int2 is the ordered interaction list (ordered according to correlation)
  # pci is bp_prior output file
  # ab is the names of the interacting bodies, can get it from the rownames and colnames of your adjacency matrix
  
  shortest.dist = rep(NA,nrow(int2))
  
  for(i in 1:nrow(int2)){
    #print(i)
    my.shortest.dist1 = my.shortest.dist2 = 0
    a = ab[int2[i,1]]
    b = ab[int2[i,2]]
    
    # if you find a in the first column, and b in the third column, you'll return shortest distance
    mm = which(pci[,1]==a)
    possib1 = any(as.character(pci[mm,3])==b)
    #cat("possib1 is",possib1,"\n")
    
    if(!is.na(possib1) && possib1==TRUE){
      ind1 = which(pci[mm,3] == b)[1] # can have multiple matches
      #cat("ind1 has length",length(ind1),"\n")
      my.shortest.dist1 = as.numeric(pci[mm[ind1],2])
    }#end if
    
    # if you find b in the first column, and a in the third column, you'll return shortest distance
    nn = which(pci[,1]==b)
    possib2 = any(as.character(pci[nn,3])==a)
    #cat("possib2 is",possib2,"\n")
    
    if(!is.na(possib2) && possib2 ==TRUE){
      ind2 = which(pci[nn,3] == a)[1] # can have multiple matches
      #cat("ind2 has length",length(ind2),"\n")
      my.shortest.dist2 = as.numeric(pci[nn[ind2],2])
    }#end if
    # can't be both true (I found some exceptions, but the shortest distance is the same)
    
    if((my.shortest.dist1 !=0 && my.shortest.dist2 !=0 ) && (my.shortest.dist1 != my.shortest.dist2)){
      warning(paste("Same match gives different shortest path in iter",i))
      shortest.dist[i] = min(my.shortest.dist1,my.shortest.dist2)
    }else if((my.shortest.dist1 !=0 && my.shortest.dist2 !=0 ) && (my.shortest.dist1 == my.shortest.dist2)){
      shortest.dist[i] = my.shortest.dist1
    }else if((my.shortest.dist1 !=0 && my.shortest.dist2 ==0) || (my.shortest.dist2 !=0 && my.shortest.dist1 ==0)){
      shortest.dist[i] = max(my.shortest.dist1,my.shortest.dist2)
    }else{}
    
  }#end for
  cat(i,"edges found\n")
  return(shortest.dist)
  
}#end function


getMutualNeighbors = function(pci,int,ab){
  N = matrix(NA,nrow=nrow(int),ncol=1,dimnames=list(c(1:nrow(int)),"neighbors"))
  
  for(i in 1:nrow(int)){
    #print(i)
    
    a = ab[int[i,1]]
    b = ab[int[i,2]]
    
    # Get the (dist1) neighbors of a from both columns 1 and 3, and then unique them
    rownums1a = intersect(which(pci[,1]==a),which(pci[,2]==1))
    rownums2a = intersect(which(pci[,3]==a),which(pci[,2]==1))
    n.a = unique(union(pci[rownums1a,3],pci[rownums2a,1]))
    
    rownums1b = intersect(which(pci[,1]==b),which(pci[,2]==1))
    rownums2b = intersect(which(pci[,3]==b),which(pci[,2]==1))
    n.b = unique(union(pci[rownums1b,3],pci[rownums2b,1]))
    
    neibor = sort(intersect(n.a,n.b))
    if(length(neibor)>0){
      N[i,1] = paste(neibor,collapse=",")
    }else{
      N[i,1] = "NA"
    }
  }#end for
  cat("Mutual neighbors added to result\n") 
  return(N)
  
}#end function

rppaCorrelNetwork = function(MAT,QUANTILE_CUTOFF,antibodyNames){
  cor1 = cor(t(MAT))
  cor2 = cor1
  cor2[lower.tri(cor2,diag=TRUE)] = 0
  
  CUTOFF = quantile(abs(cor2[cor2!=0]),QUANTILE_CUTOFF)
  int = which(abs(cor2)>CUTOFF,arr.ind=TRUE)
  cors = round(cor2[which(abs(cor2)>CUTOFF)],6)
  int.ordered = int[order(abs(cors),decreasing=TRUE),]
  #nrow(interactions)
  ab = antibodyNames
  intMat = matrix(NA,ncol=7,nrow=nrow(int.ordered))
  colnames(intMat) = c("antibody1","gene1","antibody2","gene2","correlation","shortestDist","mutualNeighbors")
  rownames(intMat) = 1:nrow(intMat)
  intMat[,1] = ab[int.ordered[,1]]
  intMat[,2] = rownames(int.ordered)
  intMat[,3] = ab[int.ordered[,2]]
  intMat[,4] = colnames(cor2)[int.ordered[,2]]
  intMat[,5] = sort(abs(cors),decreasing=TRUE)
  intMat[,6] = getShortestDistance(pci,int.ordered,ab)
  intMat[,7] = getMutualNeighbors(pci,int.ordered,ab)
  
  return(intMat)
}#end function

rppaGenenetNetwork = function(MAT,CUTOFF,antibodyNames){
  
  tmat = t(MAT)
  
  # reliable estimation of the partial correlation matrix
  inferred.pcor = ggm.estimate.pcor(tmat)
  # Estimating optimal shrinkage intensity lambda (correlation matrix): 0.1488 
  #dim(inferred.pcor) # 165 by 165
  test.results = ggm.test.edges(inferred.pcor)
  # Estimate (local) false discovery rates (partial correlations):
  # Step 1... determine cutoff point
  # Step 2... estimate parameters of null distribution and eta0
  # Step 3... compute p-values and estimate empirical PDF/CDF
  # Step 4... compute q-values and local fdr
  # Step 5... prepare for plotting
  
  signif = test.results$prob > CUTOFF
  #sum(signif)
  sigmat = test.results[signif,]
  
  # visualization
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("Rgraphviz")
  
  #node.labels = ab
  #gr = ggm.make.graph(GeneNet.int,node.labels)
  #gr
  #A graphNEL graph with undirected edges
  #Number of Nodes = 165 
  #Number of Edges = 160 
  #show.edge.weights(gr)
  #ggm.plot.graph(gr,show.edge.labels=F,layoutType="fdp")
  ab = antibodyNames      
  sigint = sigmat[,2:3]
  
  intMat = matrix(NA,ncol=7,nrow=nrow(sigint))
  colnames(intMat) = c("antibody1","gene1","antibody2","gene2","correlation","shortestDist","mutualNeighbors")
  rownames(intMat) = 1:nrow(intMat)
  intMat[,1] = ab[sigint[,1]]
  intMat[,2] = rownames(MAT)[sigmat[,2]]
  intMat[,3] = ab[sigint[,2]]
  intMat[,4] = rownames(MAT)[sigmat[,3]]
  intMat[,5] = sigmat[,1]
  intMat[,6] = getShortestDistance(pci,sigint,ab)
  intMat[,7] = getMutualNeighbors(pci,sigint,ab)
  
  return(intMat)
}#end function

rppaParcorNetwork = function(MAT,QUANTILE_CUTOFF,antibodyNames){
  covmat = cov(t(MAT))
  invmat = solve(covmat)
  invmat2 = invmat
  
  invmat2[lower.tri(invmat2,diag=TRUE)] = 0
  CUTOFF = quantile(abs(invmat2[invmat2!=0]),QUANTILE_CUTOFF)
  int = which(abs(invmat2)>CUTOFF,arr.ind=TRUE)
  
  vals = round(invmat2[which(abs(invmat2)>CUTOFF)],3)
  int.ordered = int[order(abs(vals),decreasing=TRUE),]
  #nrow(interactions)
  ab = antibodyNames
  intMat = matrix(NA,ncol=7,nrow=nrow(int.ordered))
  colnames(intMat) = c("antibody1","gene1","antibody2","gene2","correlation","shortestDist","mutualNeighbors")
  rownames(intMat) = 1:nrow(intMat)
  intMat[,1] = ab[int.ordered[,1]]
  intMat[,2] = rownames(int.ordered)
  intMat[,3] = ab[int.ordered[,2]]
  intMat[,4] = rownames(MAT)[int.ordered[,2]]
  intMat[,5] = sort(abs(vals),decreasing=TRUE)
  intMat[,6] = getShortestDistance(pci,int.ordered,ab)
  intMat[,7] = getMutualNeighbors(pci,int.ordered,ab)
  
  return(intMat)
}#end function

getSummaryMatrix = function(CORMAT,PCORMAT,GNMAT){
  summary = matrix(NA,nrow=3,ncol=7,dimnames=list(c("cor","partialCor","GeneNet"),c("dist1","dist1%","dist2","dist2%","known","known%","total")))
  
  summary[1,] = fillSummaryRow(CORMAT)
  summary[2,] = fillSummaryRow(PCORMAT)
  summary[3,] = fillSummaryRow(GNMAT)
  
  return(summary)
}#end function

fillSummaryRow = function(MAT){
  tot = nrow(MAT)
  res = matrix(NA,1,5)
  res[1] = length(which(as.numeric(MAT[,6])==1))
  res[2] = round((res[1] / tot),4)
  res[3] = length(which(as.numeric(MAT[,6])==2))
  res[4] = round((res[3] / tot),4)
  res[5] = res[1] + res[3]
  res[6] = round((res[5] / tot),4)
  res[7] = tot
  
  return(res)
}#end function