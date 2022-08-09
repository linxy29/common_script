# Xinyi Lin, 202206
# Goal: Extract gene expression matrix from a Seurat object for scVelo
# Useage: Rscript /home/linxy29/code/R/common_script/get_pySCENIC.R ../day7_preprocess/day7iPSC_wUMAP.RDS "day7"

# !!!!This isn't finished!!!! 

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
write.csv(t(as.matrix(seuratObj@assays$RNA@counts)),file = paste0(SampleName, ".csv"))
message(Sys.time(), "\tDone.")

