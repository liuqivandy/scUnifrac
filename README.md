scUnifrac
==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#example)

<a name="introduction"/>

# Introduction

scUnifrac is to quantify cell subpopulation diversity between two single-cell transcriptome profiles, which calculates the distance, estimates the statistical significance of the distance, identifies the subpopulations that are significant different between two profiles, finds genes that mark subpopulations, and predicts cell types of subpopulations. If two single-cell RNA-seq profiles have identical cell subpopulation structures (the same subpopulations with the same proportion), the distance is zero, whereas, the distance is one if completely different subpopulations.In addition to return the distance and the pvalue, scUnifrac will generate a report including figures illustrating subpopulation structure, marker genes and predicted cell types in each subpopulation. 

<a name="installation"/>

# Installation

To generate a report , you need to install pandoc (https://github.com/jgm/pandoc/releases/tag/2.2.1) . After installation , update your path to include the directory where pandoc’s binaries are installed.

<br>

To install scUnifrac, use

	source("https://bioconductor.org/biocLite.R")
	biocLite("liuqivandy/scUnifrac")
  
<a name="example"/>

# Usage

After installing scUnifrac, use the following code to run a simple example

	library(scUnifrac)
	#load two example datasets, one is from mouse colon, the other is from mouse pancreas
	load(system.file("extdata", "colon1.Rdata", package = "scUnifrac"))
	load(system.file("extdata", "pan1.Rdata", package = "scUnifrac"))
	
	##run scUnifrac on the two example datasets
	result<-scUnifrac(data1=colon1,data2=pan1)
	result
	
	##run scUnifrac, use a reference datasets to predict cell types
	#load Mouse cell atlas [Han et al., 2018, Cell 172, 1091–1107]. The atlas is used as a reference to predict cell types
	load(system.file("extdata", "ref.expr.Rdata", package = "scUnifrac"))
	result<-scUnifrac(data1=colon1,data2=pan1,ref.expr=ref.expr)
	
	#run scUnifrac on two identical samples
	ind<-sample(c(1:ncol(colon1)), ncol(colon1)/2)
	result<-scUnifrac(data1=colon1[,ind],data2=colon1[,-ind],report=F)
	
	

