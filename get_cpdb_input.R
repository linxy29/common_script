# Xinyi Lin, 202207
# Goal: Extract gene expression matrix from a Seurat object for cellphoneDB
# Useage: Rscript /home/linxy29/code/R/common_script/get_cpdb_input.R ../day7_preprocess/day7iPSC_wUMAP.RDS "day7"

options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)


### Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(EnsDb.Hsapiens.v79))

rm(list=ls())

### Functions



### Arguments to be provided when executing script
args = commandArgs(trailingOnly=TRUE)
SeuratPath <- args[1] # The path to the seurat object. It should be save as RDS file
SampleName <- args[2] # Sample name that will be used for output files

### Example parameters
#SeuratPath <- "/storage/holab/linxy/iPSC/day7_preprocess/day7iPSC_wUMAP.RDS" # The path to the seurat object. It should be save as RDS file
#SampleName <- "day7" # Sample name that will be used for output files

message(Sys.time(), "\tLoading the Seurat data.")
seuratObj <- readRDS(SeuratPath)

if (length(args) == 3){
  EnsemblePath <- args[3]
  ensembleF = fread(EnsemblePath)
  colnames(ensembleF) = c("GENEID", "SYMBOL", "info")
} else {
  geneSymbols <- rownames(seuratObj@assays$RNA)
  ensembleF <- ensembldb::select(EnsDb.Hsapiens.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
}

## generate normalized count data
message(Sys.time(), "\tGenerating the count data.")
count_raw <- seuratObj@assays$RNA@counts[ensembleF$SYMBOL,]
count_norm<- apply(count_raw, 2, function(x) (x/sum(x))*10000)
rownames(count_norm) = ensembleF[match(rownames(count_norm), ensembleF$SYMBOL),]$GENEID
count_normDF = count_norm %>% as.data.frame() %>% 
  mutate(Gene = rownames(count_norm)) %>% 
  dplyr::select(Gene, everything())
write.table(count_norm, paste0(SampleName, "_count.txt"), sep="\t", quote=F)

## generate metadata
message(Sys.time(), "\tGenerating the meta data.")
meta_data <- cbind(colnames(count_norm), seuratObj$seurat_clusters %>% as.character())
colnames(meta_data) = c("Cell", "cell_type")
write.table(meta_data, paste0(SampleName, "_meta.txt"), sep="\t", quote=F, row.names=F)
message(Sys.time(), "\tDone.")

