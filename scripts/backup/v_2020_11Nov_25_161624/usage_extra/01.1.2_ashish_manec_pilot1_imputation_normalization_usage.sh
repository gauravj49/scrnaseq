# pwd

# ############### Imputation ########################

ipython3

# Denoising using DCA
import scanpy as sc
from dca.api import dca
from gjainPyLib import *

ipython # Python 3.6.8 (default, Jan 14 2019, 11:02:34)
# aedata = sc.read_text('input/manac/T_cellranger_filtered_manac_counts_raw_genesymbols.txt')
aedata = sc.read_10x_mtx(
    '/home/rad/users/ashish/10x_947_manec_counts/outs/filtered_feature_bc_matrix/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)                                # write a cache file for faster subsequent reading
aedata.var_names_make_unique()
aedata
# Out[30]: AnnData object with n_obs × n_vars = 6398 × 31053


output_file = 'output/manac/01_preprocessed/02_denoised'
create_dir("{0}".format(get_file_info(output_file)[0]))
sc.pp.filter_genes(aedata, min_cells=10)
sc.pp.filter_cells(aedata, min_genes=2000)
aedata.var_names_make_unique()
aedata
# AnnData object with n_obs × n_vars = 5881 × 15165 
#     obs: 'n_genes' 15165
#     var: 'n_cells' 5881

# Denoise the data on raw counts
dca(aedata)
# DCA: Successfully preprocessed 15165 genes and 5881 cells.

# Save the denoised data
aedata.to_df().to_csv("{0}_T_cellranger_filtered_manac_counts_raw_genesymbols.txt".format(get_file_info(output_file)[3]), sep='\t', header=True, index=True, index_label="CellId")

aedata.to_df().T.to_csv("{0}_cellranger_filtered_manac_counts_raw_genesymbols.txt".format(get_file_info(output_file)[3]), sep='\t', header=True, index=True, index_label="GeneSymbol")

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
outputDir <- paste0(jobdir, "/output/manac/01_preprocessed"); system(paste("mkdir -p", outputDir, sep=' '))

rawinfile  <- paste0(jobdir ,"/input/manac/cellranger_filtered_manac_counts_raw_genesymbols.txt")
dcainfile  <- paste0(outputDir,"/02_denoised_cellranger_filtered_manac_counts_raw_genesymbols.txt")

# Get the raw input matrix
rawMat    <- data.frame(read.table(rawinfile, header=TRUE, sep="\t", na.strings=c(NA, NaN, Inf), row.names=1))
dcaMat    <- data.frame(read.table(dcainfile, header=TRUE, sep="\t", na.strings=c(NA, NaN, Inf), row.names=1))

# Normalize using scTransform
# Normalized data is normData$y
normRawData <- sctransform::vst(as.matrix(rawMat), return_gene_attr = FALSE)
normDCAData <- sctransform::vst(as.matrix(dcaMat))

# save the current R session to the file "RSession.Rda"
save.session("sctransform_RSession.Rda")
# # exit R without saving data
# q("no")
# # restart R
# R
# # load a saved R session from "RSession.Rda"
# restore.session("sctransform_RSession.Rda")


# Save file
normRawOutfile <- paste0(outputDir,"/03_cellranger_filtered_manac_counts_normalizedRaw_genesymbols.txt")
write.table(normRawData$y, file = normRawOutfile, row.names = T, col.names=T, sep = '\t', quote = F)
tnormRawOutfile<- paste(outputDir,"/03_T_cellranger_filtered_manac_counts_normalizedRaw_genesymbols.txt", sep="")
write.table(t(normRawData$y), file = tnormRawOutfile, row.names = T, col.names=T, sep = '\t', quote = F)

normDCAOutfile <- paste0(outputDir,"/03_cellranger_filtered_manac_counts_normalizedDCA_genesymbols.txt")
write.table(normDCAData$y, file = normDCAOutfile, row.names = T, col.names=T, sep = '\t', quote = F)
tnormDCAOutfile<- paste(outputDir,"/03_T_cellranger_filtered_manac_counts_normalizedDCA_genesymbols.txt", sep="")
write.table(t(normDCAData$y), file = tnormDCAOutfile, row.names = T, col.names=T, sep = '\t', quote = F)



