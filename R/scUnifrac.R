doReport<-function(plotData, outputFile){
  #reportRmd <- "e:/sqh/programs/scUnifrac/inst/report/scUnifracReport.Rmd"
  reportRmd <- system.file("report/scUnifracReport.Rmd", package="scUnifrac")
  
  outputFile <- getAbsolutePath(outputFile)
  output_dir = dirname(outputFile)
  output_file = basename(outputFile)
  
  cat("Output report to:", outputFile, "\n")
  rmarkdown::render(reportRmd,
                    output_dir = output_dir,
                    output_file = output_file,
                    params = list(data = plotData, dir = output_dir))
}

scUnifracFromFile<-function(outputFile, sampleFile1, sampleName1, sampleFile2, sampleName2, refExprFile, genenum=500, ncluster=10, nDim=4, normalize=T, cache=TRUE){
  if(cache){
    plotData<-prepareReportDataFromFile(sampleFile1, sampleName1, sampleFile2, sampleName2, refExprFile, genenum, ncluster, nDim, normalize, cachePrefix=outputFile)
  }else{
    plotData<-prepareReportDataFromFile(sampleFile1, sampleName1, sampleFile2, sampleName2, refExprFile, genenum, ncluster, nDim, normalize)
  }
  
  doReport(plotData, outputFile)
}

scUnifrac<-function(outputFile, data1, sampleName1, data2, sampleName2, ref.expr, genenum=500, ncluster=10, nDim=4, normalize=T, cache=TRUE){
  if(cache){
    plotData<-prepareReportData(data1, sampleName1, data2, sampleName2, ref.expr, genenum, ncluster, nDim, normalize, cachePrefix=outputFile)
  }else{
    plotData<-prepareReportData(data1, sampleName1, data2, sampleName2, ref.expr, genenum, ncluster, nDim, normalize)
  }
  doReport(plotData, outputFile)
}
