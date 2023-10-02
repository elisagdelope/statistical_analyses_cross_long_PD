# Title: func_summarize_pathway_level
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script contains a function that summarizes expression data at pathway level (across gene members, for each sample/individual) 


# Different aggregation statistics can be computed: median (default), mean, sd, min, max, pca, mds, ttest, wilcoxon test, kolmogorov test.
# Minimum size per aggregation group is by default 10, and can be changed with the parameter minsize.
# Usage: 
# summarize_pathway_level(exprssion_matrix, gene_sets, stats_type, minsize) to get the 'type' statistic of expression across gene members of the list of gene sets.
# the result is a matrix-like table with the aggregation statistics for each group, i.e. 1 row per pathway, n columns (being n the number of samples).

summarize_pathway_level = function(exprs, gsets=NULL, type="median", minsize = 10) {
  if(is.null(gsets) && is.null(database))
  {
    stop("Either the gsets or the database parameter must be specified!")
  }
  
  exprs <- as.matrix(exprs)
  genenames = rownames(exprs)
  
  # list of pathway expression matrices
  pathmat = matrix(0, nrow=length(gsets), ncol=ncol(exprs))
  rownames(pathmat) = rep("", nrow(pathmat))
  
  count = 0
  count_0 = 0
  cat('\n',length(gsets),' pathways read.\n')
  for(j in 1:length(gsets)) {
    if(j %% 100 == 0){
      print(paste("iteration",j))
    }
    
    gset = gsets[[j]]
    
    mapid = match(gset, genenames)
    notna = which(!is.na(mapid))
    
    if(length(notna) <  minsize) {
#      print(paste(names(gsets[j]), "did not pass minsize threshold; # matching genes:", length(notna)))
      count_0 = count_0 + 1
      next
    } else {
#      print(paste(names(gsets[j]),"passed minsize threshold"))
    }
    curpathmat = exprs[mapid[notna],]
    
    meanpathvec = NULL
    if(type == "mean") {
      meanpathvec = apply(curpathmat, 2, mean)
    } else if(type=="min"){
      meanpathvec = apply(curpathmat, 2, min)
    } else if(type=="max"){
      meanpathvec = apply(curpathmat, 2, max)    
    } else if(type=="sd"){
      meanpathvec = apply(curpathmat, 2, sd)
    } else if(type=="pca"){
      rem = which(apply(curpathmat, 1, var)==0)
      curpathmatfilt = curpathmat
      if(length(rem))
        curpathmatfilt = curpathmat[-rem,]
      if(length(curpathmatfilt))
      {
        pca    <- prcomp(t(curpathmatfilt), retx=T, scale=T) # scaled pca 
        scores <- pca$x[,1]
        meanpathvec = scores
      } else {
        meanpathvec = rep(0, ncol(exprs))
      }
    } else if(type=="mds") {
      meanpathvec <- as.vector(cmdscale(dist(t(curpathmat)), k = 1))
    } else if(type=="ttest"){
      print(j)
      path_outmat <- exprs[-mapid[notna],]
      path_ttestres = sapply(1:ncol(curpathmat), function(x) {dat = t.test(curpathmat[,x], path_outmat[,x], alternative="greater"); list(dat$stat, dat$p.value)})
      path_ttest = as.numeric(path_ttestres[1,])
      #path_ttestpval = as.numeric(path_ttestres[2,])
      meanpathvec = path_ttest
    } else if(type=="wilcox"){
      print(j)
      path_outmat <- exprs[-mapid[notna],]
      path_mwres = sapply(1:ncol(curpathmat), function(x) {dat = wilcox.test(curpathmat[,x], path_outmat[,x], alternative="greater"); list(dat$stat, dat$p.value)})
      path_mw = as.numeric(path_mwres[1,])
      #path_mwpval = as.numeric(path_mwres[2,])        
      meanpathvec = path_mw
    } else if(type=="kolmogorov"){
      print(j)
      path_outmat <- exprs[-mapid[notna],]
      path_ksres = sapply(1:ncol(curpathmat), function(x) {dat = ks.test(curpathmat[,x], path_outmat[,x], alternative="greater"); list(dat$stat, dat$p.value)})       
      path_ks = as.numeric(path_ksres[1,])
      #path_kspval = as.numeric(path_ksres[2,])          
      meanpathvec = path_ks
    } else {
      meanpathvec = apply(curpathmat, 2, median)
    }       
    
    count = count + 1
    pathmat[count,] = meanpathvec
    rownames(pathmat)[count] = names(gsets)[j]
  }
  print(paste0(count, " successful aggregations over minsize"))
  print(paste0(count_0, " failed aggregations under minsize"))
  
  if (count > 1) {
    pathmat = pathmat[1:count,]
    colnames(pathmat) <- colnames(exprs)
    return(pathmat)
  } else if (count == 1) {
    pathmat = data.frame(matrix(pathmat[1,], 1), row.names = rownames(pathmat)[1])
    colnames(pathmat) <- colnames(exprs)
    return(pathmat)
  }
#  return(pathmat)
}

