doReport<-function(plotData, outputFile, outputPdf=F, htmlForPrint=F){
  if(outputPdf){
    reportRmd <- system.file("report/scUnifracReportPdf.Rmd", package="scUnifrac")
  }else{
    if(htmlForPrint){
      reportRmd <- system.file("report/scUnifracReportPrint.Rmd", package="scUnifrac")
    }else{
      reportRmd <- system.file("report/scUnifracReport.Rmd", package="scUnifrac")
    }
  }

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

scUnifracFromFile<-function(outputFile, sampleFile1, sampleName1, sampleFile2, sampleName2, refExprFile, genenum=500, ncluster=10, nDim=4, normalize=T, report=T, cache=TRUE, outputPdf=FALSE, htmlForPrint=FALSE){
  if(cache){
    plotData<-prepareReportDataFromFile(sampleFile1, sampleName1, sampleFile2, sampleName2, refExprFile, genenum, ncluster, nDim, normalize, report, cachePrefix=outputFile)
  }else{
    plotData<-prepareReportDataFromFile(sampleFile1, sampleName1, sampleFile2, sampleName2, refExprFile, genenum, ncluster, nDim, normalize, report)
  }
  
  if(report){
	  doReport(plotData, outputFile, outputPdf, htmlForPrint)
  }
  
  return(list(distance=plotData$distance,
              pvalue=plotData$pvalue))
}

scUnifrac<-function(outputFile, data1, sampleName1, data2, sampleName2, ref.expr, genenum=500, ncluster=10, nDim=4, normalize=T, report=T, cache=TRUE, outputPdf=FALSE, htmlForPrint=FALSE){
  if(cache){
    plotData<-prepareReportData(data1, sampleName1, data2, sampleName2, ref.expr, genenum, ncluster, nDim, normalize, report, cachePrefix=outputFile)
  }else{
    plotData<-prepareReportData(data1, sampleName1, data2, sampleName2, ref.expr, genenum, ncluster, nDim, normalize, report)
  }

  if(report){
    doReport(plotData, outputFile, outputPdf, htmlForPrint)
  }
  
  return(list(distance=plotData$distance,
              pvalue=plotData$pvalue))
}
