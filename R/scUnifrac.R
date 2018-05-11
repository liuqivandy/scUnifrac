doReport<-function(plotData, outputFile, pdf=F){
  if(pdf){
    reportRmd <- system.file("report/scUnifracReportPdf.Rmd", package="scUnifrac")
  }else{
    reportRmd <- system.file("report/scUnifracReport.Rmd", package="scUnifrac")
  }

  outputPrefix<-sub("[.][^.]*$", "", outputFile) 
  outputRmd<-paste0(outputPrefix, ".Rmd")
  
  file.copy(reportRmd, outputRmd, overwrite=TRUE)
  
  outputFile <- getAbsolutePath(outputFile)
  cat("Output report to:", outputFile, "\n")
  
  knitr::knit_meta(class=NULL, clean = TRUE)
  rmarkdown::render(outputRmd, params = list(data = plotData, outFile = outputFile))
}

scUnifracFromFile<-function(outputFile, sampleFile1, sampleName1, sampleFile2, sampleName2, refExprFile, genenum=500, ncluster=10, nDim=4, normalize=T, cache=TRUE, pdf=FALSE){
  if(cache){
    plotData<-prepareReportDataFromFile(sampleFile1, sampleName1, sampleFile2, sampleName2, refExprFile, genenum, ncluster, nDim, normalize, cachePrefix=outputFile)
  }else{
    plotData<-prepareReportDataFromFile(sampleFile1, sampleName1, sampleFile2, sampleName2, refExprFile, genenum, ncluster, nDim, normalize)
  }
  
  doReport(plotData, outputFile, pdf)
}

scUnifrac<-function(outputFile, data1, sampleName1, data2, sampleName2, ref.expr, genenum=500, ncluster=10, nDim=4, normalize=T, cache=TRUE, pdf=FALSE){
  if(cache){
    plotData<-prepareReportData(data1, sampleName1, data2, sampleName2, ref.expr, genenum, ncluster, nDim, normalize, cachePrefix=outputFile)
  }else{
    plotData<-prepareReportData(data1, sampleName1, data2, sampleName2, ref.expr, genenum, ncluster, nDim, normalize)
  }
  doReport(plotData, outputFile, pdf)
}
