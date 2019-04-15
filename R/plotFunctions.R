
Plotdiff_tree<-function(otu.tab, tree, diffpop_perm_max,relpop_perm_max, legendtxt){
  diffp<-Diff_population(otu.tab,tree)
  diffpop<-abs(diffp$cum[,1]-diffp$cum[,2])
  relpop<-diffp$relpop
  ##the permuated proportion difference 
  
  ##find the node with signficant proprotion difference		
  diffid<- diffpop>diffpop_perm_max
  
  ### define the color###      
  col<-rep("black",nrow(tree$edge))
  col[(diffp$cum[,1]-diffp$cum[,2]>0) & (relpop>apply(cbind(0.95,relpop_perm_max),1,min))] <-"pink"
  col[(diffp$cum[,1]-diffp$cum[,2]<0) & (relpop>apply(cbind(0.95,relpop_perm_max),1,min))]<-"lightblue"
  
  col[(diffp$cum[,1]-diffp$cum[,2]>0) & diffid]<-"red"
  col[(diffp$cum[,1]-diffp$cum[,2]<0) & diffid]<-"blue"
  edgewidth<-diffpop*5
  edgewidth[edgewidth<0.5]<-0.5
  
  plot(tree,edge.color=col,edge.width=edgewidth, label.offset=1.5,cex=min(1,40/ncol(otu.tab)))
  
  pievalue<-t((otu.tab)/apply((otu.tab),1,sum))
  
  pievalue<-pievalue/apply(pievalue,1,sum)
  tiplabels(pie=pievalue,piecol=c("red","blue"),adj=1,cex=min(0.6,0.6*20/ncol(otu.tab)))
  legend("topleft",legend=c(legendtxt$samplenames), pch=16,col=c("red","blue"))
  title(paste("Distance=",round(legendtxt$dis,2),", pval=",legendtxt$pval))
}

Plotdiff_tSNE<-function(otu.tab, data, memb, diffnodes,colcluster,legendtxt ){
  tsne_out<-Rtsne(data)
  ncluster<-ncol(otu.tab)
  nsample<-apply(otu.tab,1,sum)
  
  shape<-c(rep(16,nsample[1]),rep(8,nsample[2]))
  plot(tsne_out$Y,col=colcluster[memb],pch=shape,xlab="tSNE1",ylab="tSNE2")
  
  legend("topleft",legend=legendtxt$samplenames,pch=c(16,8))
  if(length(diffnodes) > 0){
    legend("topright",legend=diffnodes,col=colcluster[diffnodes],lty=1,bty="n")
    for (i in 1:length(diffnodes)){
      cord<-apply(tsne_out$Y[memb==diffnodes[i],],2,median)
      text(cord[1],cord[2],labels=diffnodes[i],col="black",cex=2)
    }
  }
  
  title(paste("Distance=",round(legendtxt$dis,2),", pval=",legendtxt$pval))
}

###find the markers in each different cluster##########################
##normdata is the normalized and log-transformed gene expression
##memb is the membership of each cell
## diffnode is the subpopulation that focus on
##colcluster is the color for each subpopulation 
##clustertext is the information about the subpopulation
##logFCcutoff and pvalcutoff is the threshold to select marker genes
Plotdiff_Markergene<-function(normdata, memb, diffnode,colcluster,clustertext,logFCcutoff=0.5,pvalcutoff=0.01){
  
  clustid<-unique(c(diffnode,which(table(memb)>20)))  #require the cluster at least 20 cells to compare with others### if diffnode<20,the node will be still included in the differential analysis
  cellid<-memb%in%clustid		
  group<-as.factor(as.character(memb[cellid]))
  
  selected_dataall<-normdata[,cellid]
  
  memb_selecteddata<-memb[cellid]
  
  group<-relevel(group,ref=as.character(diffnode))
  design <- model.matrix(~group)
  fit <- lmFit(selected_dataall, design)
  fit <- eBayes(fit, trend=TRUE, robust=TRUE)
  uplist<-downlist<-NULL
  for (i in 2: ncol(design)){
    result<-topTable(fit,coef=i,number=nrow(normdata))
    downgene<-result[result$logFC<(-logFCcutoff) & result$adj.P.Val<pvalcutoff,]
    
    top10down<-rownames(downgene[1:min(nrow(downgene),10),])
    downlist<-c(downlist, top10down)
    upgene<-result[result$logFC>logFCcutoff & result$adj.P.Val<pvalcutoff,]
    top10up<-rownames(upgene[1:min(nrow(upgene),10),])
    uplist<-c(uplist,top10up)
  }
  
  ###generate the map
  genelist<-unique(c(uplist,downlist))
  markerdata<-selected_dataall[rownames(selected_dataall)%in%genelist,]
  colval<-rep(diffnode,sum(memb==diffnode))
  markerdata_reorder<-markerdata[,memb_selecteddata==diffnode]
  for (bindind in 1:length(clustid)){
    if (clustid[bindind]!=diffnode){
      markerdata_reorder<-cbind(markerdata_reorder, markerdata[,memb_selecteddata==clustid[bindind]])
      colval<-c(colval,rep(clustid[bindind],sum(memb==clustid[bindind])))	
    }			 
  }
  
  scaledmatrix<-t(scale(t(markerdata_reorder)))
  scaledmatrix[scaledmatrix>3]<-3
  scaledmatrix[scaledmatrix<(-3)]<-(-3)
  
  rhc<-hclust(as.dist(1-cor(t(scaledmatrix))),"ave")
  scaledmatrix_reorder<-scaledmatrix[rhc$order,]
  rownum<-nrow(scaledmatrix)
  colnum<-ncol(scaledmatrix)
  layout(matrix(c(2,4,3,1), 2, 2, byrow = TRUE),width=c(1,6),heights=c(1,8))
  par(yaxs="i")
  par(xaxs="i")
  par(mar=c(3,0,0,6))
  
  image(1:colnum,1:rownum,t(scaledmatrix_reorder),col= colorRampPalette(c("blue","white", "red"))(100),zlim=c(-3,3),axes=F,xlab="",ylab="")
  rowcex<-min(0.8, 50 / rownum * 0.8)
  axis(4,1:rownum,labels=rownames(scaledmatrix_reorder),las=2,tick=F,cex.axis=rowcex)
  mtext(clustertext,side=1,line=1,at=colnum/2,cex=1.5)
  
  par(mar=c(2,1,2,1))
  image(1:7,1,matrix(seq(-3,3,1),ncol=1),col= colorRampPalette(c("blue","white", "red"))(100),zlim=c(-3,3),axes=F)
  box(lty="solid",col="black")
  mtext("-3",side=1,line=0.5,at=1,cex=0.8)
  mtext("3",side=1,line=0.5,at=7,cex=0.8)
  
  par(mar=c(3,0,0,0))
  plot(as.dendrogram(rhc),horiz=T,leaflab="none",axes=F,xlab="",ylim=c(0,rownum))
  par(mar=c(0,0,2,6))
  plot(NA, bty="n",axes=FALSE,xlim=c(0,colnum), ylim=c(0,5),ylab="",xlab="")
  for (i in 1:colnum){
    rect(i-1,0,i,1,col=colcluster[colval[i]],border=NA)
  }
  legend("top",legend=unique(colval),fill=colcluster[unique(colval)],ncol=10,bty="n")
}	

###predict the cell type
##testdata is the gene expression of the query cells
##cluster text is the information about the subpopulation
Plotdiff_Celltype<-function(testdata, ref.expr, clustertext){
  commongene<-intersect(rownames(testdata),rownames(ref.expr))
  ##require more than 300 genes in common to predict cell types
  if (length(commongene)>300){
    tst.match <- testdata[commongene,]
    ref.match<-ref.expr[commongene,]
    
    cors <- cor(ref.match,tst.match)
    
    cors_index <- apply(cors,2,function(x){return(order(x,decreasing=T)[1:3])})
    cors_index <- sort(unique(as.integer(cors_index)))
    scblast.result <- apply(cors,2,function(x) rownames(cors)[which.max(x)])
    if (length(cors_index)>1){
      
      cors_in = cors[cors_index,]
      
      colnames(cors_in)<-colnames(testdata)
      
      rhc<-hclust(dist(t(cors_in)),method="ward.D2")
      hc<-hclust(dist((cors_in)),method="ward.D2")
      rownum<-nrow(cors_in)
      colnum<-ncol(cors_in)
      layout(matrix(c(4,2,3,1), 2, 2, byrow = TRUE),width=c(1,6),heights=c(1,6))
      par(mar=c(3,1,1,10))
      par(xaxs="i")
      par(yaxs="i")
      image(1:colnum,1:rownum,t(cors_in[hc$order,rhc$order]),col=colorRampPalette(c("gray","white","red"))(100),axes=F) 
      rowcex<-min(0.6, 50 / rownum * 0.6)
      axis(4,1:rownum,labels=rownames(cors_in)[hc$order],las=2,tick=F,cex.axis=rowcex)
      mtext(clustertext,side=1,line=1,at=colnum/2,cex=1.5)
      par(mar=c(0,1,1,10))
      plot(as.dendrogram(rhc),axes=F,leaflab="none")
      par(mar=c(3,1,1,0))
      plot(as.dendrogram(hc),horiz=T,axes=F,leaflab="none")
      par(mar=c(2,1,2,1))
      maxval<-round(max(cors_in),digits=1)
      zval<-seq(0,maxval,0.1)
      image(1:length(zval),1,matrix(zval,ncol=1),col= colorRampPalette(c("gray","white", "red"))(100),axes=F,xlab="",ylab="")
      mtext("0",side=1,line=0.5,at=1,cex=0.8)
      mtext(maxval,side=1,line=0.5,at=length(zval),cex=0.8)
      
      box(lty="solid",col="black")
      
    } else {
      
      cors_in=matrix(cors[cors_index,],nrow=1)
      colnames(cors_in)<-colnames(testdata)	
      rownames(cors_in)<-rownames(cors)[cors_index]		
      rownum<-nrow(cors_in)
      colnum<-ncol(cors_in)
      rhc<-hclust(dist(t(cors_in)),method="ward.D2")
      layout(matrix(c(3,2,0,1), 2, 2, byrow = TRUE),width=c(1,6),heights=c(1,6))
      par(mar=c(3,1,1,10))
      par(xaxs="i")
      par(yaxs="i")
      image(1:colnum,1:rownum,t(t(cors_in[,rhc$order])),col=colorRampPalette(c("gray","white","red"))(100),axes=F,xlab="",ylab="") 
      rowcex<-min(0.6, 50 / rownum * 0.6)
      axis(4,1:rownum,labels=rownames(cors_in),las=2,tick=F,cex.axis=rowcex)
      mtext(clustertext,side=1,line=1,at=colnum/2,cex=1.5)
      par(mar=c(0,1,1,10))
      plot(as.dendrogram(rhc),axes=F,leaflab="none")
      
      par(mar=c(2,1,2,1))
      maxval<-round(max(cors_in),digits=1)
      zval<-seq(0,maxval,0.1)
      image(1:length(zval),1,matrix(zval,ncol=1),col= colorRampPalette(c("gray","white", "red"))(100),axes=F,xlab="",ylab="")
      mtext("0",side=1,line=0.5,at=1,cex=0.8)
      mtext(maxval,side=1,line=0.5,at=length(zval),cex=0.8)
      
      box(lty="solid",col="black")
    }
  }
}
