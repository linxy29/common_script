# Xinyi Lin, 202206
# Goal: Extract gene expression matrix from a Seurat object for scVelo
# Useage: Rscript /home/linxy29/Code/common_script/get_scVelo_input.R ../seuratObj/2022-08-02_MES.CR.6.1.2.5000nF.clustered.RDS "ncc"

options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)


### Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratWrappers))

rm(list=ls())

### Functions


### Arguments to be provided when executing script
Path <- commandArgs(trailingOnly=TRUE)[1] # The path to the seurat object. It should be save as RDS file
SampleName <- commandArgs(trailingOnly=TRUE)[2] # Sample name that will be used for output files


seuratObj <- readRDS(Path)
message(Sys.time(), "\tLoaded the Seurat data.")

## This is for CLI pySCENIC
write.csv(Cells(seuratObj), file = paste0(SampleName, "_cellID.csv"), row.names = FALSE)
write.csv(Embeddings(seuratObj, reduction = "umap"), file = paste0(SampleName, "_embeddings.csv"))
write.csv(Idents(seuratObj), file = paste0(SampleName, "_clusters.csv"))
message(Sys.time(), "\tDone.")
