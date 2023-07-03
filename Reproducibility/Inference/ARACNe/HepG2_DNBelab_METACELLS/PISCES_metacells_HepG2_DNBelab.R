library(viper)
library(Seurat)
library(igraph)
library(biomaRt)
library(uwot)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(PISCES)

# Input matrix: UMI counts (it can be created from a .h5ad file adapting the script convert_h5_seurat.R)
pisces.obj <- readRDS("./HepG2_DNBelab_UMI.rds")

print("Loaded data and Created PISCES object")
pisces.obj <- SCTransform(pisces.obj, variable.features.n = 5000)
pisces.obj <- RunPCA(pisces.obj, verbose = FALSE)
pisces.obj <- RunUMAP(pisces.obj, dims = 1:30, verbose = FALSE)
pisces.obj <- CorDist(pisces.obj, pca.feats = 10)

my.clust <- as.character(pisces.obj[[]]$batch)
names(my.clust) <- rownames(pisces.obj[[]])

## if your data has fewer than 10K UMIs / cell, we recommend generating metacells
## the clustering vector can be any of those generated in the previous step; here we use the PISCES clustering
gexp.dist <- pisces.obj@assays$SCT@misc$dist.mat

pisces.metacell.mats <- MetaCells(as.matrix(pisces.obj@assays$RNA@counts),dist.mat = gexp.dist,clust.vect = my.clust)

dir.create('./HepG2_DNBelab_METACELLS/pisces-clust_aracne-mats_new/')
for (mcn in names(pisces.metacell.mats)) {
  saveRDS(pisces.metacell.mats[[mcn]], 
          file = paste('./HepG2_DNBelab_METACELLS/pisces-clust_aracne-mats_new/mats_pisces-', 
                       mcn, '.rds', sep = ''))
}
