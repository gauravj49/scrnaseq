# pwd

# ############### Imputation ########################

ipython3

# Denoising using DCA
import scanpy as sc
from dca.api import dca
from gjainPyLib import *

ipython # Python 3.6.8 (default, Jan 14 2019, 11:02:34)
aedata = sc.read_text('input/manac/T_cellranger_filtered_manac_counts_raw_genesymbols.txt')
aedata.var_names_make_unique()
aedata
# Out[30]: AnnData object with n_obs × n_vars = 6398 × 30958


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
normRaw_data <- sctransform::vst(as.matrix(rawMat[complete.cases(rawMat * 0), , drop=FALSE]))
normDCA_data <- sctransform::vst(as.matrix(dcaMat))

# save the current R session to the file "RSession.Rda"
save.session("RSession.Rda")

# exit R without saving data
q("no")

# restart R
R

# load a saved R session from "RSession.Rda"
restore.session("RSession.Rda")








# 2) SCNORM

# # Remove duplicates from the raw file
# awk '!seen[$1]++' input/publicDS3/12weeks_hippo_publicDS3_dzne_counts_raw_genesymbols.txt > input/publicDS3/12weeks_hippo_publicDS3_dzne_counts_raw_genesymbols_sorted.txt
output/publicDS3/scanpy/


# Run in R
R 
library(data.table)
library(SCnorm)


# Create the output directories and files
jobdir    <- "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq"
outputDir <- paste0(jobdir, "/output/manac/01_preprocessed"); system(paste("mkdir -p", outputDir, sep=' '))

rawinfile  <- paste0(jobdir ,"/input/manac/T_cellranger_filtered_manac_counts_raw_genesymbols.txt")
dcainfile  <- paste0(outputDir,"/02_denoised_cellranger_filtered_manac_counts_raw_genesymbols.txt")

# Get the raw input matrix
rawMat    <- data.frame(read.table(rawinfile, header=TRUE, sep="\t", na.strings=c(NA, NaN, Inf), row.names=1))
dcaMat    <- data.frame(read.table(dcainfile, header=TRUE, sep="\t", na.strings=c(NA, NaN, Inf), row.names=1))
# rawMatDT    <- fread(rawinfile); rawMatDF <- as.data.frame(rawMatDT); rownames(rawMatDT)    <- names(rawMat[0])
# 		sampleTable$sampleName   <- NULL
# dcaMatDT    <- fread(dcainfile); 

# Run SCnorm
rawCond    <- rep(1, dim(rawMat)[1])

countDeptEst <- plotCountDepth(Data = rawMat, Conditions = rawCond, FilterCellProportion = .1, NCores=3)
normRawMat <- SCnorm(Data=rawMat, K=1, Conditions = rawCond)


# conditions  <- rep("1", dim(rawMat)[1])
conditions <- c()
for (c in colnames(rawMat)){
	if(grepl('CA1',c,fixed=TRUE)){
		conditions <- c(conditions, 'CA1')
	}else if(grepl('CA23',c,fixed=TRUE)){
		conditions <- c(conditions, 'CA23')
	}else if(grepl('DG',c,fixed=TRUE)){
		conditions <- c(conditions, 'DG')
	}
}
normRawMat  <- SCnorm(Data=rawMat, Conditions=conditions, K=2)

conditions <- c()
for (c in colnames(dcaMat)){
	if(grepl('CA1',c,fixed=TRUE)){
		conditions <- c(conditions, 'CA1')
	}else if(grepl('CA23',c,fixed=TRUE)){
		conditions <- c(conditions, 'CA23')
	}else if(grepl('DG',c,fixed=TRUE)){
		conditions <- c(conditions, 'DG')
	}
}
normDcaMat  <- SCnorm(Data=dcaMat, Conditions=conditions)

# Get normalized counts matrices
normRawData <- results(normRawMat, type = "NormalizedData")
normDcaData <- results(normDcaMat, type = "NormalizedData")

# Save file
normRawOutfile <- paste(outputDir,"/scnorm_12weeks_hippo_publicDS3_dzne_counts_raw_genesymbols.txt", sep="")
write.table(normRawData, file = normRawOutfile, row.names = T, col.names=T, sep = '\t', quote = F)
tnormRawOutfile<- paste(outputDir,"/T_scnorm_12weeks_hippo_publicDS3_dzne_counts_raw_genesymbols.txt", sep="")
write.table(t(normRawData), file = tnormRawOutfile, row.names = T, col.names=T, sep = '\t', quote = F)

normDcaOutfile <- paste(outputDir,"/scnorm_denoised_12weeks_hippo_publicDS3_dzne_counts_raw_genesymbols.txt", sep="")
write.table(normDcaData, file = normDcaOutfile, row.names = T, col.names=T, sep = '\t', quote = F)
tnormDcaOutfile<- paste(outputDir,"/T_scnorm_denoised_12weeks_hippo_publicDS3_dzne_counts_raw_genesymbols.txt", sep="")
write.table(t(normDcaData), file = tnormDcaOutfile, row.names = T, col.names=T, sep = '\t', quote = F)



