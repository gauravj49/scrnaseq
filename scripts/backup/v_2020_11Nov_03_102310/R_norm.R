
# Create the output directories and files
rawinfile <- "output/manec/tissue/01_preprocessing/manec_tissues_merged_except1079/counts/01_raw_manec_tissues_merged_except1079_filtered.txt"

# Get the raw input matrix
rawMat    <- data.frame(read.table(rawinfile, header=TRUE, sep="\t", na.strings=c(NA, NaN, Inf), row.names=1))

conditions <- (as.vector(unlist(read.table("conditions.txt", sep="\t")))) 

normRawMat  <- SCnorm(Data=rawMat, Conditions=conditions, K=2)
normRawData <- results(normRawMat, type = "NormalizedData")

normRawOutfile <- "output/manec/tissue/01_preprocessing/manec_tissues_merged_except1079/counts/02_normalizedRaw_manec_tissues_merged_except1079_filtered.txt"
write.table(normRawData, file = normRawOutfile, row.names = T, col.names=T, sep = '\t', quote = F)




