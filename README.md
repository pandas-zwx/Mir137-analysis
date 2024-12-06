# Mir137 RNA-Seq and scRNA-Seq Analysis Scripts

This repository contains two R scripts for RNA sequencing (RNA-Seq) data analysis, including bulk RNA-Seq and single-cell RNA-Seq (scRNA-Seq). These scripts perform various analyses, such as filtering low-quality cells, identifying marker genes, and characterizing cell clusters.

Scripts Overview

1. main_137_DiffGenes.R

	•	Purpose: This script is designed for analyzing bulk RNA-Seq data, including differential gene expression analysis.
	•	Key Features:
	•	Reads in bulk RNA-Seq count data from cutdiff analysis.
	•	Performs data normalization and quality control(FPKM2TPM.R).
	•	Conducts differential gene expression analysis using statistical methods.
	•	Outputs visualizations such as heatmaps and volcano plots.

2. Mir137_scRNA.R

	•	Purpose: This script focuses on single-cell RNA-Seq (scRNA-Seq) data analysis, with specific functionalities for data preprocessing, cell clustering, and marker identification.
	•	Key Features:
	•	Filters low-quality cells based on predefined criteria (e.g., mitochondrial gene expression, total UMIs).
	•	Normalizes and scales the scRNA-Seq data.
	•	Performs principal component analysis (PCA) for dimensionality reduction.
	•	Identifies cell clusters using clustering algorithms (e.g., Louvain or Leiden algorithms).
	•	Finds marker genes for each cell cluster.
	•	Annotates cell clusters based on marker genes and biological knowledge.

Requirements

	•	R version: ≥ 4.0.0
	•	Required Packages:
	•	Seurat
	•	dplyr
	•	ggplot2
	•	Additional packages as specified in the script headers.

How to Use

	1.	Clone the repository:

git clone https://github.com/your-username/your-repository.git
cd your-repository


	2.	Install the required R packages by running:

install.packages(c("Seurat", "dplyr", "ggplot2", "DESeq2"))


	3.	Run the scripts in R or RStudio:
	•	For bulk RNA-Seq analysis:

source("main_137_DiffGenes.R")


	•	For scRNA-Seq analysis:

source("Mir137_scRNA.R")



Outputs

	•	main_137_DiffGenes.R:
	•	Differential gene expression tables.
	•	Heatmaps and volcano plots for visualization.
	•	Mir137_scRNA.R:
	•	Filtered cell datasets.
	•	PCA plots and clustering visualizations.
	•	Marker gene tables for identified clusters.

Notes

	•	Ensure that the input data is in the correct format as specified in the scripts.
	•	Modify script parameters (e.g., filtering thresholds, clustering resolution) as needed for your specific dataset.

Feel free to use or modify this README.md file as needed for your project! If you need further assistance, let me know.
