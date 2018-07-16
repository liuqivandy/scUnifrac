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
#' @param ref.expr matrix; the data matrix of the reference cell atlas; the ref.expr is used to predict the cell types of single cell from data1 and data2 (default: NULL)
#' @param genenum integer; Number of highly variable genes to build the cell population structure (default: 500)
#' @param ncluster integer; Number of clusters to divide cells  (default: 10)
#' @param nDim integer; Number of PCA dimensions to build the cell population structure  (default: 4)
#' @param normalize logical; Indicate whether normalize data1 and data2 (default: TRUE, normalize to the total count and log2 transform)
#' @param report logical;  Indicate whether to generate a report (default: TRUE); if set to FALSE, scUnifrac will not generate the report but calculate the distance and the pvalue; set to FALSE, when users have more than two samples to compare and only want to calculate the pairwise distance
#' 
#' @return List with the following elements:
#' \item{distance}{The distance of cell population diversity between two single-cell RNA-seq datasets}
#' \item{p-value}{The statistical signficance of the distance}

#' @examples
#' #library(scUnifrac)  
#' ##load the two example datasets 
#' #load(system.file("extdata", "colon1.Rdata", package = "scUnifrac"))
#' #load(system.file("extdata", "pan1.Rdata", package = "scUnifrac"))
#' 
#' #result<-scUnifrac( data1=colon1, data2=pan1)  #this function will return distance and pvalue between data1 and data2 and will generate a report in the work directory
#' #result
#' ##load the mouse cell altas from Han et al., 2018, Cell 172, 1091–1107. The atlas is used as a reference to predict cell types of data1 and data2
#' #load(system.file("extdata", "ref.expr.Rdata", package = "scUnifrac"))
#' #result<-scUnifrac( data1=colon1, data2=pan1,ref.expr=ref.expr, outputFile="scUnifrac_report.html") #the report includes the predicted cell types of each cell in each differential subpopulation between data1 and data2
#' 
#' ##test two samples with similar populations
#' #ind<-sample(c(1:ncol(colon1)), ncol(colon1)/2)
#' #result<-scUnifrac(data1=colon1[,ind], data2=colon1[,-ind],ref.expr=ref.expr, outputFile="scUnifrac_report.html")
#' 
#' @import limma ape permute GUniFrac Rtsne R.utils knitr kableExtra rmdformats statmod
#' @importFrom devtools session_info
#' 
#' @export

scUnifrac<-function(data1, sampleName1="S1", data2, sampleName2="S2", ref.expr=NULL, genenum=500, ncluster=10, nDim=4, normalize=T, report=T, outputFile="scUnifrac_report.html"){
  plotData<-prepareReportData(data1, sampleName1, data2, sampleName2, ref.expr, genenum, ncluster, nDim, normalize, report)

  if(report){
    doReport(plotData, outputFile)
  }
  
  return(list(distance=plotData$distance,
              pvalue=plotData$pvalue))
}


#' scUnifrac_multi 
#' 
#' @description Quantify cell population diversity among single cell RNA-seq (>= two datasets)


#' @param dataall matrix; the combined data matrix of all datasets, row is the gene symbol, the column is the cell id 
#' @param group a vector; giving the sample id/name to which each cell belongs to
#' @param genenum integer; Number of highly variable genes to build the cell population structure (default: 500)
#' @param ncluster integer; Number of clusters to divide cells  (default: 10)
#' @param nDim integer; Number of PCA dimensions to build the cell population structure  (default: 4)
#' @param normalize logical; Indicate whether normalize data1 and data2 (default: TRUE, normalize to the total count and log2 transform)

#' @return List with the following elements:
#' \item{distance}{The pairwise distance matrix of cell population diversity among single-cell RNA-seq datasets}
#' \item{p-value}{The statistical signficance matrix of the distance}

#' @examples
#' #library(scUnifrac)  
#' ##load the two example datasets 
#' #load(system.file("extdata", "colon1.Rdata", package = "scUnifrac"))
#' #load(system.file("extdata", "pan1.Rdata", package = "scUnifrac"))
#' ##split the colon data into two datasets
#' #colon1_1<-colon1[,1:500]
#' #colon1_2<-colon1[,501:1000]
#' #result<-scUnifrac_multi(dataall=cbind(colon1_1,colon1_2,pan1),group=c(rep("c1",500),rep("c2",500),rep("pan",ncol(pan1))))  #this function will return distance and pvalue between data1 and data2 and will generate a report in the work directory
#' #result

#' @import limma ape permute GUniFrac Rtsne R.utils knitr kableExtra rmdformats statmod
#' @importFrom devtools session_info
#' 
#' @export

scUnifrac_multi<-function(dataall,group,genenum=500,ncluster=10,nDim=4,normalize=T){
    
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
    
    hc<-hclust(dist(hvgdata_pca),"ave")
    memb <- cutree(hc, k = ncluster)
       
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
    return(list(dis=dist.obs, counttable=count.table, pvalue=pvalue))
}


