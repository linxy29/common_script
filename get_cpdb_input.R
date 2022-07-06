# Xinyi Lin, 202207
# Goal: Extract gene expression matrix from a Seurat object for pySCENIC
# Useage: Rscript /home/linxy29/code/R/common_script/get_pySCENIC.R ../day7_preprocess/day7iPSC_wUMAP.RDS "day7"

# !!!!This isn't finished!!!! 

options(max.print = 500)
options(stringsAsFactors = FALSE)
options(scipen = 999)


### Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

rm(list=ls())

### Functions



### Arguments to be provided when executing script
SeuratPath <- commandArgs(trailingOnly=TRUE)[1] # The path to the seurat object. It should be save as RDS file
EnsemblePath <- commandArgs(trailingOnly=TRUE)[2] # The path to the features file
SampleName <- commandArgs(trailingOnly=TRUE)[3] # Sample name that will be used for output files

### Example parameters
SeuratPath <- commandArgs(trailingOnly=TRUE)[1] # The path to the seurat object. It should be save as RDS file
EnsemblePath <- commandArgs(trailingOnly=TRUE)[2] # The path to the features file
SampleName <- commandArgs(trailingOnly=TRUE)[3] # Sample name that will be used for output files

seuratObj <- readRDS(Path)
message(Sys.time(), "\tLoaded the Seurat data.")

ensembleF = fread(str_c(pathSPC, "/GCTB/cellranger_output/L56_v3/outs/filtered_feature_bc_matrix/features.tsv.gz"), header = FALSE)
colnames(ensembleF) = c("ensembleID", "geneID", "info")

load(str_c(pathSPC, "/GCTB/GCTB6/GCTB6.combined.rdata"))

## generate normalized count data
GCTB6_row <- GCTB6.combined@assays$RNA@counts
GCTB6_norm<- apply(GCTB6_row, 2, function(x) (x/sum(x))*10000)
rownames(GCTB6_norm) = ensembleF[match(rownames(GCTB6_norm), ensembleF$geneID),]$ensembleID
GCTB6_normDF = GCTB6_norm %>% as.data.frame() %>% 
  mutate(Gene = rownames(GCTB6_norm)) %>% 
  dplyr::select(Gene, everything())
#write.table(GCTB6_normDF[,1:101],"D:/Data/GCTB/cpDB/GCTB6_count_test2.txt", row.names=FALSE, sep="\t", quote=F)
write.table(GCTB6_norm[,1:100],str_c(pathSPC, "/GCTB/cpDB/GCTB6_count_test.txt"), sep="\t", quote=F)
write.table(GCTB6_norm,str_c(pathSPC, "/GCTB/cpDB/GCTB6_count.txt"), sep="\t", quote=F)

## generate metadata
GCTB6_meta_data <- cbind(colnames(GCTB6_norm), GCTB6.combined$seurat_clusters %>% as.character())
colnames(GCTB6_meta_data) = c("Cell", "cell_type")
GCTB6_meta_data_test <- cbind(colnames(GCTB6_norm[,1:100]), GCTB6.combined$seurat_clusters %>% as.character() %>% head(100) %>% str_c("cluster", .))
colnames(GCTB6_meta_data_test) = c("Cell", "cell_type")
#write.table(GCTB6_meta_data_test , "D:/Data/GCTB/cpDB/GCTB6_meta_test2.txt", sep="\t", quote=F, row.names=F, col.names = TRUE)
#write.table(GCTB6_meta_data_test , "D:/Data/GCTB/cpDB/GCTB6_meta_test.txt", sep="\t", quote=F, row.names=F)
write.table(GCTB6_meta_data, str_c(pathSPC, "/GCTB/cpDB/GCTB6_meta.txt"), sep="\t", quote=F, row.names=F)
message(Sys.time(), "\tDone.")

