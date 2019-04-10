
Diff_population<-function(otu.tab, tree){
  
  otu.tab <- as.matrix(otu.tab)
  row.sum <- rowSums(otu.tab)
  otu.tab <- otu.tab / row.sum
  n <- nrow(otu.tab)
  tip.label <- tree$tip.label
  otu.tab <- otu.tab[, tip.label]
  
  ntip <- length(tip.label)
  nbr <- nrow(tree$edge)	
  edge <- tree$edge
  edge2 <- edge[, 2]
  br.len <- tree$edge.length
  
  #  Accumulate OTU proportions up the tree	
  cum <- matrix(0, nbr, n)							# Branch abundance matrix
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i)
    cum[tip.loc, ] <- cum[tip.loc, ] + otu.tab[, i]	
    node <- edge[tip.loc, 1]						# Assume the direction of edge 
    node.loc <- which(edge2 == node)
    while (length(node.loc)) {
      cum[node.loc, ] <- cum[node.loc, ] + otu.tab[, i]		
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node)
    }
  }
  relpop=cum[,1]/apply(cum,1,sum)
  relpop=ifelse(relpop>0.5,relpop,1-relpop)
  return(list(cum=cum,relpop=relpop))
}

.getCacheFile<-function(cachePrefix){
  return(paste0(cachePrefix, "_scUnifrac.rdata"))
}

prepareReportData<-function(data1, sampleName1="S1", data2, sampleName2="S2", cluster=NULL, ref.expr=NULL, genenum=500, ncluster=10, nDim=4, normalize=T, report=T, cachePrefix){  
  if(missing(ref.expr) ){
    ref.expr<-NULL
  }
  
  saveCache<-!missing(cachePrefix)
  
  if(saveCache){
    cacheFile = .getCacheFile(cachePrefix)
    if(file.exists(cacheFile)){
      cat("Loading from cache file", cacheFile, " ...\n")
      load(cacheFile)
      return(plotData)
    } 
  }
  
  samplenames=c(sampleName1,sampleName2)
    
  if (is.null(colnames(data1))) {colnames(data1)<-1:ncol(data1)}
  if (is.null(colnames(data2))) {colnames(data2)<-1:ncol(data2)}
  if (is.null(samplenames)) {samplenames<-c("S1","S2")}
  
  nobs1<-ncol(data1)
  nobs2<-ncol(data2)
  
  data_all<-cbind(data1,data2)
  if (normalize){
    sumdata<-apply(data_all,2,sum)
    #normalize the data    
    normdata<-t(log2(t(data_all)/sumdata*median(sumdata)+1))
  } else {
    normdata<-data_all
  }
  #find the highly variable genes	
  var_data<-apply(normdata,1,var)
  mean_data<-apply(normdata,1,mean)
  
  reorderid<-order(var_data,decreasing=T)
  
  hvgdata<-normdata[reorderid[1:genenum],]
  
  ##perform PCA, keep the nDim   
  pcaresult<-prcomp(t(hvgdata))
  hvgdata_pca<-pcaresult$x[,1:nDim]
  
  if (is.null(cluster)){
  ##build the hierichical cluster 
       hc<-hclust(dist(hvgdata_pca),"ave")
       memb <- cutree(hc, k = ncluster)
	  } else { if (length(cluster)!=(nobs1+nobs2)) 
		      stop ('the length of cluster information should be equal to the cell numbers')
		   memb<-as.integer(as.factor(cluster))
		   ncluster<-length(unique(cluster))
	  }
	  
  
  ##generate the table 
  count.table<-matrix(0,nrow=2,ncol=ncluster)
  colnames(count.table)<-1:ncluster
  gr1count<- table(memb[1:nobs1])
  gr2count<-table(memb[(nobs1+1):(nobs1+nobs2)])
  
  count.table[1,names(gr1count)]<-gr1count
  count.table[2,names(gr2count)]<-gr2count
  
  cent <- NULL
  for(k in 1:ncluster){
    cent <- rbind(cent, colMeans(hvgdata_pca[memb == k, , drop = FALSE]))
  }
  
  hc1 <- hclust(dist(cent), method = "ave", members = table(memb))
  tree1<-as.phylo(hc1)
  
  ##calculate the distance  
  unifracs <- GUniFrac(count.table, tree1, alpha=c(0, 0.5, 1))$unifracs
  dist.obs<-unifracs[, , "d_1"][1,2]
  
  ##permuate the samples, calculate the pvalue   
  nperm=1000
  control <- how(nperm = nperm, within = Within(type = "free"))
  t.permu <- numeric(length = control$nperm) +1
  diffpop_perm<-relpop_perm<-NULL
  for(numPer in seq_along(t.permu)) {
    want <- permute(numPer, nobs1+nobs2, control)
    ## calculate distance
    count.table.perm<-matrix(0,nrow=2,ncol=ncluster)
    colnames(count.table.perm)<-1:ncluster
    gr1count.perm<- table(memb[want[1:nobs1]])
    gr2count.perm<-table(memb[want[(nobs1+1):(nobs1+nobs2)]])
    count.table.perm[1,names(gr1count.perm)]<-gr1count.perm
    count.table.perm[2,names(gr2count.perm)]<-gr2count.perm
    
    unifracs.perm <- GUniFrac(count.table.perm, tree1, alpha=c(0, 0.5, 1))$unifracs
    
    diffp_perm<-Diff_population(count.table.perm, tree1)
    diffpop_perm<-cbind(diffpop_perm,abs(diffp_perm$cum[,1]-diffp_perm$cum[,2]))
    relpop_perm<-cbind(relpop_perm,diffp_perm$relpop)
    t.permu[numPer]<- unifracs.perm[, , "d_1"][1,2]
  }
  
  pvalue<-sum(dist.obs<t.permu)/nperm
  
  if(!report){
	return(list(distance=dist.obs,
	            pvalue=pvalue))
  }
  
  diffpop_perm_max<-apply(diffpop_perm,1,max)
  relpop_perm_max<-apply(relpop_perm,1,max)
  
  ##legned txt for the following tree and tSNE figures
  legendtxt<-list(samplenames=samplenames,dis=dist.obs,pval=pvalue)
  
  ##find the node with signficant proprotion difference		
  diffp<-Diff_population(count.table, tree1)
  diffpop<-abs(diffp$cum[,1]-diffp$cum[,2])
  relpop<-diffp$relpop
  diffind<- diffpop>diffpop_perm_max

  ###plot the subpopulation difference on tSNE plot
  diffnodes<-tree1$edge[diffind,2]
  diffnodes<-diffnodes[diffnodes<=ncluster]
  #diffnodes sorted by the difference
  diffval<-sort(abs(count.table[1,diffnodes]/nobs1-count.table[2,diffnodes]/nobs2),decreasing=T)
  diffnodes<-as.numeric(names(diffval))
  
  ####generate the color for each cluster with different subpopulation
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  colcluster=sample(color,ncluster)
  notdiffnodes<-!rep(1:ncluster)%in% diffnodes
  colcluster[notdiffnodes]<-"gray95"
  
  plotData<-list(distance=dist.obs,
	               pvalue=pvalue,
				         ref.expr=ref.expr,
                 count.table=count.table, 
                 normdata=normdata,
                 tree=tree1, 
                 diffind=diffind, 
                 diffpop_perm_max=diffpop_perm_max,
                 relpop_perm_max=relpop_perm_max, 
                 hvgdata_pca=hvgdata_pca, 
                 memb=memb,
                 diffnodes= diffnodes,
                 diffval=diffval,
                 colcluster=colcluster,
                 samplenames=samplenames,
                 legendtxt=legendtxt)
  
  if(saveCache){
    save(plotData, file=cacheFile)
  }
  
  return(plotData)
}
