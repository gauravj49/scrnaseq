library(reticulate)
#py_install("anndata")
library(sceasy)
library(Seurat)

h5ad_file <- "05_SCSA_Leiden_Clustering_tregCNS_adata.h5ad"
sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",outFile='05_SCSA_Leiden_Clustering_tregCNS_adata.rds')

# X -> counts
# An object of class Seurat
# 15835 features across 20149 samples within 1 assay
# Active assay: RNA (15835 features, 0 variable features)
# 3 dimensional reductions calculated: pca, tsne, umap
scsa_seurat <- readRDS("05_SCSA_Leiden_Clustering_tregCNS_adata.rds")

# > scsa_seurat
# An object of class Seurat
# 15835 features across 20149 samples within 1 assay
# Active assay: RNA (15835 features, 0 variable features)
# 3 dimensional reductions calculated: pca, tsne, umap
