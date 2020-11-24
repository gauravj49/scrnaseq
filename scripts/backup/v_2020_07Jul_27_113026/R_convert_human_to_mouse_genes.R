# SOURCE: https://rjbioinformatics.com/2016/10/14/converting-mouse-to-human-gene-names-with-biomart-package/

require("biomaRt")
library(data.table)

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
 
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
 
humanx <- unique(genesV2[, 2])
 
# Print the first 6 genes found to the screen
print(head(humanx))
return(humanx)
}

hgDT     <- fread("regev_lab_cell_cycle_genes_hsa.txt", header=F)
humGenes <- hgDT$V1

mouseGenes <- convertHumanGeneList(humGenes)
