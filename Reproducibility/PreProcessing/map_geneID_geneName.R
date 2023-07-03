library(Seurat)
library(tidyverse)
library(hablar)
library(data.table)
library(biomaRt)
source("/Users/jonathan/Desktop/IIT/Cerase_single_cell/ANALYSIS/Seurat_analysis/myfunctions.R")

ensembl107 <- useEnsembl(biomart = 'genes', 
                         dataset = 'hsapiens_gene_ensembl',
                         version = 107)

mapGeneIds <- function(counts){
  my.attributes <- c("ensembl_gene_id","external_gene_name")
  my.query <- getBM(attributes = my.attributes,
                    values = rownames(counts), 
                    mart = ensembl107)
  my.query<- my.query[my.query$ensembl_gene_id %in% row.names(counts),]
  counts <- counts[rownames(counts) %in% my.query$ensembl_gene_id,]
  my.query <- my.query[match(rownames(counts), my.query$ensembl_gene_id),]
  print(length(rownames(counts)))
  print(length(my.query$ensembl_gene_id))
  all(my.query$ensembl_gene_id==rownames(counts))
  rownames(counts) <- make.unique(my.query$external_gene_name)
  counts <- counts[rownames(counts)!='',]
  print(head(rownames(counts)))
  counts
}


options(future.globals.maxSize = 1000 * 1024^2)
wd <- "/Users/jonathan/Desktop/IIT/INTERACTomics/scRNA-seq_data/K562/"
setwd(wd)

# Load counts for SMART-seq3
counts.smart3 <- read.csv("./Smart-seq3/K562_Smart_seq3_umi_counts.txt",sep='\t',row.names = 1, header = TRUE)
counts.smart3 <- mapGeneIds(counts.smart3)
write.csv("./Smart-seq3/K562_Smart_seq3_umi_counts_gnames.txt")