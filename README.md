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

You need to install some bioconductor required libraries first:

	source("https://bioconductor.org/biocLite.R")
	biocLite("limma")
	
Then you can install scUnifrac by:

	library(devtools)
	devtools::install_github("liuqivandy/scUnifrac")
  
<a name="example"/>

# Quick start

Here we show the most basic steps.

	library(scUnifrac)
	load(system.file("extdata", "colon1.Rdata", package = "scUnifrac"))
	load(system.file("extdata", "pan1.Rdata", package = "scUnifrac"))
	load(system.file("extdata", "ref.expr.Rdata", package = "scUnifrac"))
	scUnifrac("scUnifracReport.pdf", colon1, "Colon", pan1, "Pan", ref.expr, cache=TRUE, pdf=TRUE )
	scUnifrac("scUnifracReport.html", colon1, "Colon", pan1, "Pan", ref.expr, cache=TRUE, pdf=FALSE )

The colon1 and pan1 dataset are two gene expression data matrix in which rownames are gene symbols and columns are samples. The rownames of two matrix should be identical since the data matrix will be merged together for analysis.
