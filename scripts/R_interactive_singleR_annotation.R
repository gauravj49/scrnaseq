# From scanpy
adata.to_df().to_csv("/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/tregCNS/individual_samples/S509/data/04_markerGenes_S509_adata.matrix", index=True, header=True, sep="\t")
# [1402 rows x 17606 columns]

adata.to_df().T.to_csv("/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/tregCNS/individual_samples/S509/data/04_markerGenes_S509_adata_T.matrix", index=True, header=True, sep="\t")
# [17606 rows x 1402 columns]


# BiocManager::install("SingleR")
suppressMessages(require(SingleR))
suppressMessages(require(SummarizedExperiment))
suppressMessages(require(data.table))

#  ==================================== ========== ======================== ========= ======================================== ==================== ==================== 
#           Retrieval function           Organism   Cell type focus           Samples                Sample types                No. of main labels   No. of fine labels  
#  ==================================== ========== ======================== ========= ======================================== ==================== ==================== 
#   HumanPrimaryCellAtlasData()          human      Non-specific                 713   microarrays of sorted cell populations                   37                  157  
#   BlueprintEncodeData()                human      Non-specific                 259   RNA-seq                                                  24                   43  
#   DatabaseImmuneCellExpressionData()   human      Immune                      1561   RNA-seq                                                   5                   15  
#   NovershternHematopoieticData()       human      Hematopoietic & Immune       211   microarrays of sorted cell populations                   17                   38  
#   MonacoImmuneData()                   human      Immune                       114   RNA-seq                                                  11                   29  
#   ImmGenData()                         mouse      Hematopoietic & Immune       830   microarrays of sorted cell populations                   20                  253  
#   MouseRNAseqData()                    mouse      Non-specific                 358   RNA-seq                                                  18                   28  
#  ==================================== ========== ======================== ========= ======================================== ==================== ==================== 

# Load singleR reference datasets for mouse
ImmGenRef      <- ImmGenData()
MouseRNAseqRef <- MouseRNAseqData()

# Read the count matrix file
inputMatrixFile <- "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/tregCNS/individual_samples/S509/data/04_markerGenes_S509_adata_T.matrix"
adataMatDT <- fread(inputMatrixFile, header=TRUE, sep="\t")
adataMatDF <- as.data.frame(adataMatDT)
rownames(adataMatDF) <- adataMatDF$V1; adataMatDF$V1 <- NULL
# Source: https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# NOTE: The keyword assay is important to explicitely mention the counts
adataSE <- SummarizedExperiment(assay=as.matrix(adataMatDF))

# Annotate it with SingleR
MouseRNAseqRefadataSE <- SingleR(test = adataSE, ref = MouseRNAseqRef, assay.type.test=1, labels = MouseRNAseqRef$label.main)
ImmGenRefadataSE      <- SingleR(test = adataSE, ref = ImmGenRef, assay.type.test=1, labels = ImmGenRef$label.main)


# # To check:
# > table(MouseRNAseqRefadataSE$labels)

#        Astrocytes           B cells   Dendritic cells Endothelial cells
#                15                 1                45                23
#  Epithelial cells       Fibroblasts      Granulocytes       Macrophages
#                84                 7                13                67
#         Microglia         Monocytes          NK cells  Oligodendrocytes
#                14               755                12               234
#           T cells
#               132


# > table(ImmGenRefadataSE$labels)

#           B cells                DC Endothelial cells  Epithelial cells
#                 1                41                25                61
#       Fibroblasts               ILC       Macrophages        Mast cells
#                26                89               683                 1
#         Microglia         Monocytes       Neutrophils          NK cells
#                79               313                13                 3
#               NKT     Stromal cells           T cells               Tgd
#                24                15                23                 5

# Save output as text files
outputImmGenAnnFile   <- "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/tregCNS/individual_samples/S509/data/05_markerGenes_S509_SingleR_ImmGenRef.txt"
outputMouseRnaAnnFile <- "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/tregCNS/individual_samples/S509/data/05_markerGenes_S509_SingleR_MouseRNAseqRef.txt"

write.table(data.frame("cellIDs" = rownames(ImmGenRefadataSE),ImmGenRefadataSE)          , file = outputImmGenAnnFile  , row.names = FALSE, col.names=T, sep = '\t', quote = F)
write.table(data.frame("cellIDs" = rownames(MouseRNAseqRefadataSE),MouseRNAseqRefadataSE), file = outputMouseRnaAnnFile, row.names = FALSE, col.names=T, sep = '\t', quote = F)
