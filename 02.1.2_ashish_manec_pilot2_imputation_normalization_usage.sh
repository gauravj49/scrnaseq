# pwd

# inputdir="/media/rad/SSD1/temp_manec"
# ############### Imputation ########################

ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)

# Denoising using DCA
import scanpy as sc
from dca.api import dca
from gjainPyLib import *

projName        = "bulk997_mouse"
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed"; create_dir("{0}".format(output_dir))
minGenesPerCell = 200
minCellsPergene = 3

# aedata = sc.read_text('input/manec/T_cellranger_filtered_manec_counts_raw_genesymbols.txt')
adata = sc.read_10x_mtx(
    '/media/rad/SSD1/temp_manec/bulk997_mouse/outs/filtered_feature_bc_matrix/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)                                # write a cache file for faster subsequent reading
adata.var_names_make_unique()
adata
# View of AnnData object with n_obs × n_vars = 4277 × 31053 , var: 'gene_ids', 'feature_types'

# Save the raw data as tab separated matrix 
adata.to_df().to_csv("{0}/01_raw_T_{1}_cellranger_filtered_manec_counts_genesymbols.txt".format(output_dir, projName), sep='\t', header=True, index=True, index_label="CellId")
adata.to_df().T.to_csv("{0}/01_raw_{1}_cellranger_filtered_manec_counts_genesymbols.txt".format(output_dir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol")


# Perform basic prefiltering
sc.pp.filter_cells(adata, min_genes=minGenesPerCell) 
sc.pp.filter_genes(adata, min_cells=minCellsPergene)
adata.var_names_make_unique()
adata
# AnnData object with n_obs × n_vars = 5881 × 15165 
#     obs: 'n_genes' 15165
#     var: 'n_cells' 5881

# Denoise the data on raw counts
dca(adata)
# DCA: Successfully preprocessed 15165 genes and 5881 cells.

# Save the denoised data
adata.to_df().to_csv("{0}/02_denoised_DCA_T_{1}_cellranger_filtered_manec_counts_genesymbols.txt".format(output_dir, projName), sep='\t', header=True, index=True, index_label="CellId")
adata.to_df().T.to_csv("{0}/02_denoised_DCA_{1}_cellranger_filtered_manec_counts_genesymbols.txt".format(output_dir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol")

# clit + D + D
#

# ############### Normalization ########################

# 1) scTransform

R 

library(data.table)
library(sctransform)
library(session)

# Create the output directories and files
jobdir    <- "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq"
outputDir <- paste0(jobdir, "/output/manec/pilot2/01_preprocessed"); system(paste("mkdir -p", outputDir, sep=' '))
projName  <- "bulk997_mouse"

rawinfile  <- paste0("/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/01_raw_bulk997_mouse_cellranger_filtered_manec_counts_genesymbols.txt")
dcainfile  <- paste0("/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/02_denoised_DCA_bulk997_mouse_cellranger_filtered_manec_counts_genesymbols.txt")

# Get the raw input matrix
rawMat    <- data.frame(read.table(rawinfile, header=TRUE, sep="\t", na.strings=c(NA, NaN, Inf), row.names=1))
dcaMat    <- data.frame(read.table(dcainfile, header=TRUE, sep="\t", na.strings=c(NA, NaN, Inf), row.names=1))

# Normalize using scTransform
# Normalized data is normData$y
normRawData <- sctransform::vst(as.matrix(rawMat), return_gene_attr = FALSE)
normDCAData <- sctransform::vst(as.matrix(dcaMat))

# save the current R session to the file "RSession.Rda"
sessionfile  <- paste0("/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/sctransform_bulk997_mouse_RSession.Rda")
save.session(sessionfile)
# # exit R without saving data
# q("no")
# # restart R
# R
# # load a saved R session from "RSession.Rda"
# restore.session("sctransform_RSession.Rda")


# Save file
normRawOutfile <- paste0(outputDir,"/03_normalizedRawCounts_", projName, "_cellranger_filtered_manec_counts_genesymbols.txt")
write.table(normRawData$y, file = normRawOutfile, row.names = T, col.names=T, sep = '\t', quote = F)
tnormRawOutfile<- paste(outputDir,"/03_normalizedRawCounts_T_", projName, "_cellranger_filtered_manec_counts_genesymbols.txt", sep="")
write.table(t(normRawData$y), file = tnormRawOutfile, row.names = T, col.names=T, sep = '\t', quote = F)

normDCAOutfile <- paste0(outputDir,"/03_normalizedDCACounts_", projName, "_cellranger_filtered_manec_counts_genesymbols.txt")
write.table(normDCAData$y, file = normDCAOutfile, row.names = T, col.names=T, sep = '\t', quote = F)
tnormDCAOutfile<- paste(outputDir,"/03_normalizedDCACounts_T_", projName, "_cellranger_filtered_manec_counts_genesymbols.txt", sep="")
write.table(t(normDCAData$y), file = tnormDCAOutfile, row.names = T, col.names=T, sep = '\t', quote = F)



