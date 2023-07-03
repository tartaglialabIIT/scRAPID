# This script takes as input the .h5ad files containing the processed scRNA-seq data 
# from Scanpy (see the PreProcessing folder) and converts them first to Seurat objects
# and next to rds objects

# Set your working directory
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(Seurat)
library(SeuratDisk)

datasets <- c("HepG2_Smartseq2_RBP_RNA500processed_HepG2_Smartseq2_RBP_RNA500",
              "HepG2_9CL_SCAN_RBP_RNA500processed_HepG2_9CL_SCAN_RBP_RNA500",
              "K562_CELseq_RBP_RNA500processed_K562_CELseq_RBP_RNA500",
              "K562_STORMseq_RBP_RNA500processed_K562_STORMseq_RBP_RNA500",
              "K562_Smartseq3_RBP_RNA500processed_K562_Smartseq3_RBP_RNA500",
              "K562_9CL_SCAN_RBP_RNA500processed_K562_9CL_SCAN_RBP_RNA500",
              "K562_UMI200_SCAN_RBP_RNA500processed_K562_UMI200_SCAN_RBP_RNA500")
identities <- c("HepG2_Smartseq2",
                "HepG2_9CL_SCAN",
                "K562_CELseq",
                "K562_STORMseq",
                "K562_Smartseq3",
                "K562_9CL_SCAN",
                "K562_UMI200_SCAN")

length(datasets)
length(identities)

dir.create("./converted_objects/")
dir.create("./rds/")

i <- 1 
for(d in datasets){
  f.name0 <- paste("./converted_objects/",d,sep = "")
  f.name1 <- paste("./processed_data/",d,sep = "")
  f.name1 <- paste(f.name1,".h5ad",sep = "")
  f.name2 <- paste(f.name0,".h5Seurat",sep = "")
  f.name3 <- paste("./rds/",d,sep = "")
  f.name3 <- paste(f.name3,".rds",sep = "")
  Convert(f.name1, dest = "h5seurat", overwrite = TRUE)
  seuratObject <- LoadH5Seurat(f.name2)
  Idents(object = seuratObject) <- "batch"
  seuratObject <- subset(x = seuratObject, idents = identities[i])
  saveRDS(seuratObject,f.name3)
  i <- i+1
}