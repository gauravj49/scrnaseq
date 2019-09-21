# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

# Generate the counts matrix from 10x genomics cell ranger and copy to input

# Define the output file
OFILE="input/manec/filtered_feature_bc_matrix.txt"

# Get the tab separated file with headers
ipython3
input_file = "input/manec/filtered_feature_bc_matrix.txt"
TinputDF   = pd.read_csv(input_file, comment='#', delimiter=",", index_col=0, engine='python')
TinputDF.to_csv('input/manec/filtered_feature_bc_matrix.txt', sep='\t', header=True, index=True, index_label='GeneID')
cltr+d+d

# ADD gene Symbols to the file
R

# load libraries
suppressPackageStartupMessages(library("biomaRt", warn.conflicts=FALSE, quietly=TRUE))
inputfile <- "input/manec/filtered_feature_bc_matrix.txt"
outputDir <- "input/manec"

# Get the input dataframe
sdata <- read.table(inputfile, header=TRUE, sep="\t")

# # Get the common gene names
# The "listDatasets" function will give you the list of all the species available (mart datasets) for a given mart:
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
head(listAttributes(ensembl))
geneSymbols <- getBM(filters= "ensembl_gene_id", attributes=c('ensembl_gene_id', 'external_gene_name'), mart = ensembl, values=levels(sdata$GeneID))

# # https://www.r-bloggers.com/converting-gene-names-in-r-with-annotationdbi/
# # remember to install it if you don't have it already
# library("AnnotationDbi")
# library("org.Mm.eg.db") 
# geneSymbols <- mapIds(org.Mm.eg.db, keys = levels(sdata$GeneID), keytype = "ENSEMBL", column="SYMBOL")

# Convert to Dataframe
genesDF           <- as.data.frame(geneSymbols)
names(genesDF)    <- c('GeneID','GeneSymbol')

# Add Gene symbol to the final df and save it
# Merge the two data frames
fdf           <- merge(genesDF, sdata, "GeneID")
rownames(fdf) <- fdf$GeneID
fdf$GeneID    <- NULL

## Writing results to table
soutputFile <- paste(outputDir,"/cellranger_filtered_manec_counts_raw_genesymbols.txt", sep="")
write.table(x=fdf, file=soutputFile, sep="\t", row.names=FALSE, quote=FALSE)
# quit()

# Transpose the raw genesymbols file for genes
datamash -H transpose < input/manec/cellranger_filtered_manec_counts_raw_genesymbols.txt > input/manec/T_cellranger_filtered_manec_counts_raw_genesymbols.txt