
##### TCGA-LUAD Data Download #####
# install.packages("tidyverse")
install.packages("BiocManager")
BiocManager::install("remotes")
BiocManager::install("ExperimentHub")
BiocManager::install(c("SparseArray", "DelayedArray", "SummarizedExperiment"), force = TRUE)
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")


setwd("TCGAdata")
library(tidyverse)# Load the tidyverse package for data manipulation and visualization
library(BiocManager) # Load the BiocManager package to manage Bioconductor packages
library(TCGAbiolinks) # Load the TCGAbiolinks package for interacting with TCGA data

#TCGA abbr of cancers：https://www.jianshu.com/p/3c0f74e85825
# Specify the cancer type, in this case, Lung Adenocarcinoma (LUAD)
cancer_type = "TCGA-LUAD"  
# A function to create a query for the Genomic Data Commons (GDC).
expquery <- GDCquery(project = cancer_type,
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts"
)
GDCdownload(expquery,directory = "GDCdata")
# GDCprepare: Prepares the downloaded data for analysis.
# summarizedExperiment: If TRUE, the data is returned as a SummarizedExperiment object.

expquery2 <- GDCprepare(expquery, directory = "GDCdata", summarizedExperiment = T)
save(expquery2, file = "luad.gdc_2022.rda") 

gene_id <- expquery2@rowRanges@elementMetadata@listData[["gene_id"]]
gene_name <- expquery2@rowRanges@elementMetadata@listData[["gene_name"]]
gene_type <- expquery2@rowRanges@elementMetadata@listData[["gene_type"]]
gene_annotation <- data.frame(gene_id,gene_name,gene_type)