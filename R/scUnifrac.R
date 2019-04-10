doReport<-function(plotData, outputFile){
  reportRmd <- system.file("report/scUnifracReportPrint.Rmd", package="scUnifrac")
  reportContentRmd <- system.file("report/scUnifracReportContent.Rmd", package="scUnifrac")

  outputPrefix<-sub("[.][^.]*$", "", outputFile) 
  outputRmd<-paste0(outputPrefix, ".Rmd")
  outputContentRmd<-paste0(outputFile, "_content.Rmd")
  
  file.copy(reportRmd, outputRmd, overwrite=TRUE)
  file.copy(reportContentRmd, outputContentRmd, overwrite=TRUE)

  outputFile <- getAbsolutePath(outputFile)
  cat("Output report to:", outputFile, "\n")
  
  knitr::knit_meta(class=NULL, clean = TRUE)
  rmarkdown::render(outputRmd, params = list(data = plotData, outFile = outputFile, contentFile=outputContentRmd))
  
  if(file.exists(outputFile)){
    file.remove(outputRmd)
    file.remove(outputContentRmd)
  }
}

#' scUnifrac 
#' 
#' @description Quantify cell population diversity between two single cell RNA-seq datasets

#' @param outputFile a character string giving the name of the report;(default: "scUnifrac_report.html") 
#' @param data1 matrix; the data matrix of the first dataset, row is the gene symbol, the column is the cell id 
#' @param sampleName1 a character string giving the name of the first dataset; (default: "S1")
#' @param data2 matrix; the data matrix of the second dataset, row is the gene symbol, the column is the cell id; data1 and data1 should have the same gene symbols
#' @param sampleName2 a character string giving the name of the second dataset; (default: "S2")
#' @param cluster a vector of cluster ID of each cell in both data1 and data2; (default: NULL); If NULL, hierarchical clustering will be used to assign cluster ID. 
#' @param ref.expr matrix; the data matrix of the reference cell atlas; the ref.expr is used to predict the cell types of single cell from data1 and data2 (default: NULL)
#' @param genenum integer; Number of highly variable genes to build the cell population structure (default: 500)
#' @param ncluster integer; Number of clusters to divide cells  (default: 10)
#' @param nDim integer; Number of PCA dimensions to build the cell population structure  (default: 4)
#' @param normalize logical; Indicate whether normalize data1 and data2 (default: TRUE, normalize to the total count and log2 transform)
#' @param report logical;  Indicate whether to generate a report (default: TRUE); if set to FALSE, scUnifrac will not generate the report but calculate the distance and the pvalue; set to FALSE, when users have more than two samples to compare and only want to calculate the pairwise distance
#' 
#' @return List with the following elements:
#' \item{distance}{The distance of cell population diversity between two single-cell RNA-seq datasets}
#' \item{pvalue}{The statistical signficance of the distance}

#' @examples
#' library(scUnifrac)  
#' ##load the two example datasets 
#' load(system.file("extdata", "colon1.Rdata", package = "scUnifrac"))
#' load(system.file("extdata", "pan1.Rdata", package = "scUnifrac"))

#' ##Calculate the distance and pvalue between data1 and data2 and generate a report summarizing the result in the work directory
#' result<-scUnifrac( data1=colon1, data2=pan1)  
#' result

#' ##load the mouse cell altas from Han et al., 2018, Cell 172, 1091–1107. The atlas is used as a reference to predict cell types of data1 and data2
#' load(system.file("extdata", "ref.expr.Rdata", package = "scUnifrac"))

#' ##Generate a report which also includes the predicted cell types by mapping each cell to the reference 
#' result<-scUnifrac( data1=colon1, data2=pan1,ref.expr=ref.expr, outputFile="scUnifrac_report.html") 
#' 
#' ##test two samples with similar populations
#' ind<-sample(c(1:ncol(colon1)), ncol(colon1)/2)
#' result<-scUnifrac(data1=colon1[,ind], data2=colon1[,-ind],ref.expr=ref.expr, outputFile="scUnifrac_report.html")
#' 
#' @import limma ape permute GUniFrac Rtsne R.utils knitr kableExtra rmdformats statmod
#' @importFrom devtools session_info
#' 
#' @export

scUnifrac<-function(data1, sampleName1="S1", data2, sampleName2="S2", cluster=NULL, ref.expr=NULL, genenum=500, ncluster=10, nDim=4, normalize=T, report=T, outputFile="scUnifrac_report.html"){
  plotData<-prepareReportData(data1, sampleName1, data2, sampleName2, cluster, ref.expr, genenum, ncluster, nDim, normalize, report)

  if(report){
    doReport(plotData, outputFile)
  }
  
  return(list(distance=plotData$distance,
              pvalue=plotData$pvalue))
}


#' scUnifrac_multi 
#' 
#' @description Quantify pairwise cell population diversity between multiple (>=2) single cell RNA-seq datasets


#' @param dataall matrix; the combined data matrix of all datasets, row is the gene symbol, the column is the cell id 
#' @param group a vector; giving the sample id/name to which each cell belongs to
#' @param cluster a vector of cluster ID of each cell belongs to; (default: NULL); If NULL, hierarchical clustering will be used to assign cluster ID.
#' @param genenum integer; Number of highly variable genes to build the cell population structure (default: 500)
#' @param ncluster integer; Number of clusters to divide cells  (default: 10)
#' @param nDim integer; Number of PCA dimensions to build the cell population structure  (default: 4)
#' @param normalize logical; Indicate whether normalize data1 and data2 (default: TRUE, normalize to the total count and log2 transform)

#' @return List with the following elements:
#' \item{distance}{The pairwise distance matrix of cell population diversity among single-cell RNA-seq datasets}
#' \item{pvalue}{The statistical signficance matrix of the distance}


#' @examples
#' library(scUnifrac)  
#' ##load the two example datasets 
#' load(system.file("extdata", "colon1.Rdata", package = "scUnifrac"))
#' load(system.file("extdata", "pan1.Rdata", package = "scUnifrac"))
#' ##generate two datasets from the colon data
#' colon1_1<-colon1[,1:500]
#' colon1_2<-colon1[,501:1000]
#' ## run scUnifrac_multi on three datasets, two from the colon, one from the pancreas
#' result<-scUnifrac_multi(dataall=cbind(colon1_1,colon1_2,pan1),group=c(rep("c1",500),rep("c2",500),rep("pan",ncol(pan1))))  
#' result

#' @import limma ape permute GUniFrac Rtsne R.utils knitr kableExtra rmdformats statmod
#' @importFrom devtools session_info
#' 
#' @export

scUnifrac_multi<-function(dataall,group,cluster=NULL,genenum=500,ncluster=10,nDim=4,normalize=T){
    
    ncell<-ncol(dataall)
    if (is.null(colnames(dataall))) {colnames(dataall)<-1:ncell}
    if (is.null(group) ) {stop('group information should be provided')}
    if (length(group)!=ncell){stop('the length of group does not match the number of cells')}
    nsample<-length(unique(group))
    if (nsample<2) {stop('should have greater than or equal to two single cell samples')}

    if (normalize){
        sumdata<-apply(dataall,2,sum)
        #normalize the data
        normdata<-t(log2(t(dataall)/sumdata*median(sumdata)+1))
    } else normdata<-dataall
    #find the highly variable genes
    var_data<-apply(normdata,1,var)
    mean_data<-apply(normdata,1,mean)
    reorderid<-order(var_data,decreasing=T)
    hvgdata<-normdata[reorderid[1:genenum],]
    
    ##perform PCA, keep the nDim
    pcaresult<-prcomp(t(hvgdata))
    hvgdata_pca<-pcaresult$x[,1:nDim]
    ##build the hierichical cluster
    if (is.null(cluster)){
          hc<-hclust(dist(hvgdata_pca),"ave")
          memb <- cutree(hc, k = ncluster)
        } else { if (length(cluster)!=dim(data)[2]) 
               stop ('the length of cluster' should be equal to the cell numbers')
                 memb<-as.integer(as.factor(cluster))
                 ncluster<-length(unique(cluster))
                }
        
       
    ##generate the table
    count.table<-table(group,memb)
    cent <- NULL
     
    for(k in 1:ncluster){
        cent <- rbind(cent, colMeans(hvgdata_pca[memb == k, , drop = FALSE]))
    }
    hc1 <- hclust(dist(cent), method = "ave", members = table(memb))
    tree1<-as.phylo(hc1)
    ##calculate the distance
    unifracs <- GUniFrac(count.table, tree1, alpha=c(0, 0.5, 1))$unifracs
    dist.obs<-unifracs[, , "d_1"]
    rownames(dist.obs)<-colnames(dist.obs)<-rownames(count.table)
    
    ##permuate the samples, calculate the pvalue
    nperm=1000
    control <- how(nperm = nperm, within = Within(type = "free"))
    t.permu <- numeric(length = control$nperm) +1
    dis.permu.array<-array(0,dim=c(nsample,nsample,nperm))
    for(numPer in seq_along(t.permu)) {
        
        want <- permute(numPer, ncell, control)
        group.permu<-group[want]
        
        ## calculate distance
        count.table.perm<-table(group.permu,memb)
        colnames(count.table.perm)<-1:ncluster
        unifracs.perm <- GUniFrac(count.table.perm, tree1, alpha=c(0, 0.5, 1))$unifracs
        dis.permu.array[,,numPer]<- unifracs.perm[, , "d_1"]
    }
    
    pvalue<-matrix(apply(apply(dis.permu.array,3,">=",dist.obs),1,sum),nrow=nsample,byrow=F)/nperm
    rownames(pvalue)<-colnames(pvalue)<-rownames(count.table)

    ##### function output####################
    return(list(distance=dist.obs, counttable=count.table, pvalue=pvalue,phylo=tree1))
}

#' scUnifrac_predictCelltype 
#' 
#' @description blast query cells against reference datasets to predict cell types


#' @param qdata matrix; the gene-cell matrix of query cells, row is the gene symbol, the column is the cell id
#' @param ref.expr matrix; the data matrix of the reference cell atlas;
#' @param anntext character; the annotation of query cells (default: "Query")
#' @param normalize logical; If TRUE, qdata will be normalized and log2transformed; (default: FALSE);if qdata is raw count table, set normalize=TRUE; 
#' @param corcutoff double; the cutoff of correlation values to predict cell types; if one reference cell type has correlation values > corcutoff with query cells and ranked the top 3, report this cell type; (default: 0; report the top 3 most correlated cell types without the requirement of correlation values)
#' @return the correlation values between query cells and the most correlated reference cell types with (cor>corcutoff):


#' @examples
#' library(scUnifrac)  
#' load one dataset 
#' load(system.file("extdata", "colon1.Rdata", package = "scUnifrac"))
#' ##load the mouse cell altas from Han et al., 2018, Cell 172, 1091–1107. 
#' load(system.file("extdata", "ref.expr.Rdata", package = "scUnifrac"))

## run scUnifrac_predictCelltype for the 100 cells of colon1 dataset
#' scUnifrac_predictCelltype(qdata=colon1[,1:100], ref.expr=ref.expr,anntext="colon",normalize=T) 
#' 
#' 
#' @export

scUnifrac_predictCelltype<-function(qdata, ref.expr, anntext="Query",normalize=FALSE, corcutoff=0){
    if (normalize){
       sumqdata<-apply(qdata,2,sum)
      #normalize the data    
       normdata<-t(log2(t(qdata)/sumqdata*10000+1))     
     }
    commongene<-intersect(rownames(qdata),rownames(ref.expr))
  ##require more than 300 genes in common to predict cell types
  if (length(commongene)>300){
    tst.match <- qdata[commongene,]
    ref.match<-ref.expr[commongene,]
    
    cors <- cor(ref.match,tst.match)
    
    cors_index <- unlist(apply(cors,2,function(x){cutoffind<-tail(sort(x),3)>corcutoff;return(order(x,decreasing=T)[1:3][cutoffind])}))
    #cors_index <- (apply(cors,2,function(x){return(order(x,decreasing=T)[1:3])}))
    cors_index <- sort(unique(as.integer(cors_index)))
    #scblast.result <- apply(cors,2,function(x) rownames(cors)[which.max(x)])
    if (length(cors_index)>1){
      
      cors_in = cors[cors_index,]
      
      colnames(cors_in)<-colnames(qdata)
      
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
      mtext(anntext,side=1,line=1,at=colnum/2,cex=1.5)
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
      
    } else if (length(cors_index)==1) {
      
      cors_in=matrix(cors[cors_index,],nrow=1)
      colnames(cors_in)<-colnames(qdata)	
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
      mtext(anntext,side=1,line=1,at=colnum/2,cex=1.5)
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
   return(cors_in)
}
