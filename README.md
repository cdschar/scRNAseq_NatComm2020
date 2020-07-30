scRNAseq\_NatComm2020
================
Chris Scharer and Tian Mi

# Repository Info

Rscripts used for processing scRNA-seq and bulk RNA-seq data and making
figures associated with the following publication:

Scharer & Patterson et al. Antibody-secreting cell destiny emerges
during the initial stages of B cell activation. Nature Communications
2020

Citation:

# single cell RNA-seq Scripts

  - scRNAseq.R: Main scRNAseq analysis pipeline based on Monocle2.
    Inlcude importing matrix from 10X CellRanger output, QC of cell
    detection, data normalization, dpFeature feature selection, t-SNE
    dimension reduction, density peaks clustering and single cell
    trajectory analysis.

  - magic\_transform.R: Script to do MAGIC imputation on both WT and
    IRF4 KO data sets.

  - scRNAseq\_lib.R: Customized plotting functions based on Monocle2.
    Functions inlcude t-SNE plots with different color schemes, t-SNE
    plots using MAGIC transformed data, SCT plot using different color
    schemas, expression pattern plot over pseudotime.

  - scRNAseq\_replotting: Usage examples of data loading and various
    plotting functions using the above library.

  - KNN\_predict.R: Script to perform KNN prediction based on bulk RNA
    seq data.

  - running\_scenic.R: Scripts to run SCENIC on scRNA data, including
    functions to plot heat maps and perform differential tests between
    score results for various clusters identified in Monocle2.

  - plot.GSEA.pathway.R: script to replot the GSEA output for figures.

### Software Versions

  - R 3.4.0
  - Monocle 2.9.0
  - cellrangerRkit 2.0.0
  - dplyr 0.7.7
  - ggplot2 3.3.2
  - viridis 0.5.1
  - tibble 1.4.2
  - Rmagic 2.0.3
  - FNN 1.1.3
  - preprocessCore 1.40.0
  - SCENIC 1.1.1-10
  - Biobase 2.38.0
  - AUCell 1.7.1

# bulk RNA-seq Scripts

  - RNAseq.sample.manifest.txt: key file for all samples and contains
    metadata for each

  - RNAseq.pipeline.R: script for initial data organization, fastq
    qc/trimming, mapping, duplicate
marking

## coverage scripts extract gene coverage, ercc coverage, and normalize the data

  - geneCounts.exon.PE.R: extracts coverage for all Entrez gene exons
    and normalizes data

  - ERCC\_coverage.R: extracts coverage for the ERCC spike in controls

  - geneTsNorm.updated.R: normalizes FPKM data to mRNAs/cells based on
    ERCC coverage

## Differential analysis, PCA, and plotting scripts

  - diff.glm\_v3.R: pairwise differential analysis with edgeR

  - basicPlots.R: function for making barplots

  - makeBarPlots.R: driver script for making barplots and grouping

  - makeMPCBarPlots.R: driver script to plot the total mRNAs for each
    sample from “Mouse.total.mpc.txt” file.

  - Mouse.total.mpc.txt: companion file to above script

  - pca.confidenceCircles.R: plot PCA using differentia

## Common functions written for many data analysis routines

  - bisTools.R: common functions used for processing files

### Software Versions

  - R 3.5.2
  - skewer version: 0.2.2
  - FastQC v0.11.4
  - Tophat v2.0.13
  - samtools 1.9
  - PICARD v1.127
  - GenomicRanges v1.34
  - edgeR v3.24.3
  - vegan v2.5.5
