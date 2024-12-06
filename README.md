# Mir137 RNA-Seq, scRNA-Seq, ATAC-Seq, and ChIP-Seq Analysis Scripts and Data

This repository contains R scripts and analysis outputs for various sequencing data, including bulk RNA-Seq, single-cell RNA-Seq (scRNA-Seq), ATAC-Seq, and ChIP-Seq. These resources are aimed at performing data preprocessing, quality control, feature annotation, motif discovery, and clustering analysis.

---

## Scripts Overview

### 1. **main_137_DiffGenes.R**
#### Purpose
This script is designed for analyzing bulk RNA-Seq data, including differential gene expression analysis.

#### Key Features
- Reads in bulk RNA-Seq count data from cutdiff.
- Performs data normalization and quality control.
- Conducts differential gene expression analysis using statistical methods.
- Outputs visualizations such as heatmaps and volcano plots.

---

### 2. **Mir137_scRNA.R**
#### Purpose
This script focuses on single-cell RNA-Seq (scRNA-Seq) data analysis, with specific functionalities for data preprocessing, cell clustering, and marker identification.

#### Key Features
- Filters low-quality cells based on predefined criteria (e.g., mitochondrial gene expression, total UMIs).
- Normalizes and scales the scRNA-Seq data.
- Performs principal component analysis (PCA) for dimensionality reduction.
- Finds marker genes for each cell cluster.
- Annotates cell clusters based on marker genes and biological knowledge.

---

## Data Files Overview

### 1. **cKO_specific_peaks_annotated.xls**
#### Description
This file contains the annotation of ATAC peaks specific to **Mir137 cKO** samples, generated using **HOMER**. It provides information about genomic regions and features associated with these peaks.

### 2. **homerResults.html**
#### Description
This file contains results from a **de novo motif search** performed using **HOMER**. It identifies novel motifs enriched in the provided ATAC-Seq data.

### 3. **knownResults.html**
#### Description
This file contains results from a **known motif search** performed using **HOMER**. It identifies motifs from a known database that are enriched in the provided ATAC-Seq data.

### 4. **Pu.1_target_broad_annotated.xls**
#### Description
This file provides the annotation of **PU.1 peaks** identified in ChIP-Seq data. It highlights genomic regions and features associated with PU.1 binding sites.

---

## Requirements

- **R version**: ≥ 4.0.0
- **Required Packages**:
  - `Seurat`
  - `dplyr`
  - `ggplot2`
  - `HOMER` (for motif discovery and ATAC-Seq/ChIP-Seq annotation)
  - Additional packages as specified in the script headers.

---

## How to Use

### 1. Clone the repository
To start, clone the repository to your local machine:
```bash
git clone https://github.com/your-username/your-repository.git
cd your-repository

### 2. Install required R packages

Make sure the necessary R packages are installed before running the scripts:
```R
install.packages(c("Seurat", "dplyr", "ggplot2"))

### 3. Run the scripts in R or RStudio

Follow the steps below to execute the scripts:
	•For bulk RNA-Seq analysis:
```R
source("main_137_DiffGenes.R")


	•For scRNA-Seq analysis:
```R
source("Mir137_scRNA.R")

### Outputs

RNA-Seq

	•Differential gene expression tables.
	•Heatmaps and volcano plots.

scRNA-Seq

	•Filtered cell datasets.
	•PCA plots and clustering visualizations.
	•Marker gene tables for identified clusters.

ATAC-Seq and ChIP-Seq

	•Annotated ATAC peaks for Mir137 cKO.
	•De novo and known motif discovery results from HOMER.
	•Annotated PU.1 peaks from ChIP-Seq.

Notes

	•Ensure that the input data is in the correct format as specified in the scripts and HOMER documentation.
	•Modify script parameters (e.g., filtering thresholds, clustering resolution) as needed for your specific dataset.
