    library(ape)
    library(permute)
    library(GUniFrac)
	library(limma)
	library(Rtsne)


	
##load three datasets, ref.expr is the mouse reference cell atlas, colon1 and pan1 are two example datasets

## Mouse cell atlas reference
load("ref.expr.Rdata")
# Two example datasets
load("colon1.Rdata")
load("pan1.Rdata")

## 	

##row is the genename, column is the cell id
##data1 and data2 should have the same rownames , cannot be NA 
##genenum is the number of highly variable genes to perform the PCA 
#ncluster is the number of subpopulations
#nDim is the number of PCA dimensions kept to identify the subpopulation
#samplenames is the names of two scRNA-seq samples
#normalize=T, then normalize to the total counts and logtransform
#plotpopulation=T, plot the gene signatures and predict cell type of each differential population
scUnifrac<-function(data1,data2,genenum=500,ncluster=10,nDim=4,normalize=T, samplenames=NULL, plotpopulation=F){
	
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
   } else normdata<-data_all
#find the highly variable genes	
    var_data<-apply(normdata,1,var)
    mean_data<-apply(normdata,1,mean)

    reorderid<-order(var_data,decreasing=T)
  
    hvgdata<-normdata[reorderid[1:genenum],]
	
##perform PCA, keep the nDim   
    pcaresult<-prcomp(t(hvgdata))
	hvgdata_pca<-pcaresult$x[,1:nDim]
		

   
   
   
##build the hierichical cluster 

    hc<-hclust(dist(hvgdata_pca),"ave")
    memb <- cutree(hc, k = ncluster)
	 

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
	
   
        diffpop_perm_max<-apply(diffpop_perm,1,max)
		relpop_perm_max<-apply(relpop_perm,1,max)
        
        ##legned txt for the following tree and tSNE figures
        legendtxt<-list(samplenames=samplenames,dis=dist.obs,pval=pvalue)
 
###plot the subpopulation difference on the tree
        pdf(paste(paste(samplenames,collapse="&"),"tree.pdf",sep="_"))
         diffind<-Plotdiff_tree(count.table,tree1,diffpop_perm_max,relpop_perm_max,legendtxt)
		 dev.off()

#####		
		
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
				 
##################################
		pdf(paste(paste(samplenames,collapse="&"),"tSNE.pdf",sep="_"))
		Plotdiff_tSNE(hvgdata_pca, memb, count.table, diffnodes,colcluster,legendtxt )
        dev.off()
	if (plotpopulation){
###find the marker genes in each cluster with subpopulation difference
		for (i in 1:length(diffnodes)){
		pdf(paste("markergene_",diffnodes[i],".pdf",sep=""))
		clustertext<-paste(paste("Subpopulation",diffnodes[i],collapse=" "),paste(count.table[,diffnodes[i]],collapse=":"),sep=",")
		Plotdiff_Markergene(normdata, memb, diffnodes[i],colcluster,clustertext)
        dev.off()
###predict the cell type of each cluster based on mouse cell atlas ##need input the unnormalized data
	    pdf(paste("celltype_",diffnodes[i],".pdf",sep=""))
		Plotdiff_Celltype(normdata[,memb==diffnodes[i]],clustertext)
        dev.off()
     }
	}
##### function output####################
	return(list(dis=dist.obs, counttable=count.table, pvalue=pvalue, diffval=diffval))
}
 


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

 	plot(tree,edge.color=col,edge.width=edgewidth)
        
    pievalue<-t((otu.tab)/apply((otu.tab),1,sum))
		
    pievalue<-pievalue/apply(pievalue,1,sum)
    tiplabels(pie=pievalue,piecol=c("red","blue"),offset=0.8,cex=0.6)
    legend("topleft",legend=c(legendtxt$samplenames), pch=16,col=c("red","blue"))
    title(paste("Distance=",round(legendtxt$dis,2),", pval=",legendtxt$pval))
	return(diffid)
		
}

Plotdiff_tSNE<-function(data, memb, otu.tab, diffnodes,colcluster,legendtxt ){

        tsne_out<-Rtsne(data)
        ncluster<-ncol(otu.tab)
        nsample<-apply(otu.tab,1,sum)
                   
        shape<-c(rep(16,nsample[1]),rep(8,nsample[2]))
        plot(tsne_out$Y,col=colcluster[memb],pch=shape,xlab="tSNE1",ylab="tSNE2")
        
        legend("topright",legend=diffnodes,col=colcluster[diffnodes],lty=1,bty="n")
        legend("topleft",legend=legendtxt$samplenames,pch=c(16,8))
        title(paste("Distance=",round(legendtxt$dis,2),", pval=",legendtxt$pval))
        
        for (i in 1:length(diffnodes)){
		   cord<-apply(tsne_out$Y[memb==diffnodes[i],],2,median)
		   text(cord[1],cord[2],labels=diffnodes[i],col="black",cex=2)
		   }
	
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
		
		###
	
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
		axis(4,1:rownum,labels=rownames(scaledmatrix_reorder),las=2,tick=F,cex.axis=0.8)
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
Plotdiff_Celltype<-function(testdata,clustertext){

	
    commongene<-intersect(rownames(testdata),rownames(ref.expr))
	##require more than 300 genes in common to predict cell types
	if (length(commongene)>300){
		matchid<-match(commongene,rownames(testdata))
	
   
		tst.match <- testdata[matchid[!is.na(matchid)],]
		matchid<-match(commongene,rownames(ref.expr))
	
		ref.match<-ref.expr[matchid[!is.na(matchid)],]
   
 
   
  
		cors <- cor(ref.match,tst.match)

		cors_index <- apply(cors,2,function(x){return(order(x,decreasing=T)[1])})
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
				axis(4,1:rownum,labels=rownames(cors_in),las=2,tick=F,cex.axis=0.6)
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
				axis(4,1:rownum,labels=rownames(cors_in),las=2,tick=F,cex.axis=0.6)
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


###



