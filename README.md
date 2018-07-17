scUnifrac
==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#example)
* [Report example](https://rawgit.com/wiki/liuqivandy/scUnifrac/scUnifrac_colon_pan.print.html)

<a name="introduction"/>

# Introduction

scUnifrac is a R package to quantify cell population diversity between two single-cell transcriptome profiles, which calculates the distance, estimates the statistical significance of the distance, identifies the subpopulations that are significant different between two profiles, finds genes that mark subpopulations, and predicts cell types of subpopulations. If two single-cell RNA-seq profiles have identical cell subpopulation structures (the same subpopulations with the same proportion), the distance is zero, whereas, the distance is one if completely different subpopulations.
<br>
Besides calculating the distance and the pvalue, scUnifrac will generate a [report](https://rawgit.com/wiki/liuqivandy/scUnifrac/scUnifrac_colon_pan.print.html) , which include results summary and figures illustrating subpopulation structure, marker genes and predicted cell types in each subpopulation. 

<a name="installation"/>

# Installation

To generate a report , you need to install pandoc (https://github.com/jgm/pandoc/releases/tag/2.2.1) . After installation , update your path to include the directory where pandoc’s binaries are installed.

<br>

To install scUnifrac, use

	source("https://bioconductor.org/biocLite.R")
	biocLite("liuqivandy/scUnifrac")
  
<a name="example"/>

# Usage

After installing scUnifrac, use the following code to run examples

	library(scUnifrac)
	#load two example datasets, one is from mouse colon, the other is from mouse pancreas
	load(system.file("extdata", "colon1.Rdata", package = "scUnifrac"))
	load(system.file("extdata", "pan1.Rdata", package = "scUnifrac"))
	
	##run scUnifrac on the two example datasets without generating a report
	result<-scUnifrac(data1=colon1,data2=pan1,report=F)
	result
	
	##run scUnifrac, use a reference datasets to predict cell types and generate a report in the work directory
	#load Mouse cell atlas [Han et al., 2018, Cell 172, 1091–1107]. The atlas is used as a reference to predict cell types
	load(system.file("extdata", "ref.expr.Rdata", package = "scUnifrac"))
	result<-scUnifrac(data1=colon1,data2=pan1,ref.expr=ref.expr,report=T)
	
	#run scUnifrac on two similar samples
	ind<-sample(c(1:ncol(colon1)), ncol(colon1)/2)
	result<-scUnifrac(data1=colon1[,ind],data2=colon1[,-ind],report=F)
	
	#run scUnifrac on multiple (>=2) scRNA-seq samples
	#generate a simulated dataset including three samples
	combineddata<-cbind(colon1[,1:500],pan1[,1:500],colon1[,501:1000])
	group<-c(rep("colon1_1",500),rep("pan1",500),rep("colon1_2",500))
	result<-scUnifrac_multi(dataall=combineddata,group=group)
	
