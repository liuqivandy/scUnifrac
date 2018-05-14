scUnifrac
==========
* [Introduction](#introduction)
* [Download and installation](#download)
* [Quick Start](#example)

<a name="introduction"/>

# Introduction

scUnifrac is to quantify cell subpopulation diversity between two single-cell transcriptome profiles, which calculates the distance, estimates the statistical significance of the distance, identifies the subpopulations that are significant different between two profiles, finds genes that mark subpopulations, and predicts cell types of subpopulations. If two single-cell RNA-seq profiles have identical cell subpopulation structures (the same subpopulations with the same proportion), the distance is zero, whereas, the distance is one if completely different subpopulations.

<a name="download"/>

# Download and installation

[optional] In order to generate report in pdf format, you need to [install MikTex and pandoc](http://rprogramming.net/create-html-or-pdf-files-with-r-knitr-miktex-and-pandoc/) first. After installation of MikTex and pandoc, we recommend you to restart your computer before you test the scUnifrac package.

<br>

Then you need to install bioconductor library limma:

	source("https://bioconductor.org/biocLite.R")
	biocLite("limma")
	
Finally you can install scUnifrac by:

	library(devtools)
	devtools::install_github("liuqivandy/scUnifrac")
  
<a name="example"/>

# Quick start

Here we show the most basic steps.

	library(scUnifrac)
	
	load(system.file("extdata", "colon1.Rdata", package = "scUnifrac"))
	load(system.file("extdata", "pan1.Rdata", package = "scUnifrac"))
	load(system.file("extdata", "ref.expr.Rdata", package = "scUnifrac"))
	
	#test two different profiles
	scUnifrac("scUnifrac_colon_pan.pdf", colon1, "Colon", pan1, "Pan", ref.expr, cache=TRUE, outputPdf=TRUE )
	scUnifrac("scUnifrac_colon_pan.html", colon1, "Colon", pan1, "Pan", ref.expr, cache=TRUE, outputPdf=FALSE )
	
	#test two similar profiles
	s1<-sample(c(1:ncol(colon1)), ncol(colon1)/2)
	s1<-s1[order(s1)]
	s2<-c(1:ncol(colon1))[!(c(1:ncol(colon1)) %in% s1)]

	data1<-colon1[,s1]
	sampleName1<-"Colon1"
	data2<-colon1[,s2]
	sampleName2<-"Colon2"
	scUnifrac("scUnifrac_colon_colon.pdf", colon1, "Colon", pan1, "Pan", ref.expr, cache=TRUE, outputPdf=TRUE )
	scUnifrac("scUnifrac_colon_colon.html", colon1, "Colon", pan1, "Pan", ref.expr, cache=TRUE, outputPdf=FALSE )

The colon1 and pan1 dataset are two gene expression data matrix in which rownames are gene symbols and columns are samples. The rownames of two matrix should be identical since the data matrix will be merged together for analysis.
