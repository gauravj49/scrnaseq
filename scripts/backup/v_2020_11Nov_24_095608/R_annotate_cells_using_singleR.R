#! /usr/bin/env Rscript

# ****************************************************************
# GOAL  : Annotate cells with singleR 
# USAGE : Rscript R_annotate_cells_using_singleR.R <input_matrix_file> <output_annotation_file>
# ****************************************************************

################ load main library ###############
suppressPackageStartupMessages(library("argparse"))

################ Main logic ###################
main <- function(){    

  ################ Load libraries ###############
  load_libraries()

  ##### Parse command line arguments ############
  args        <- check_options()
  inputfile   <- args$inputfile
  outputfile  <- args$outputfile

  cat("- Input command line arguments ...\n")
  cat(paste0("\t- inputfile   = ",inputfile,"\n"))
  cat(paste0("\t- outputfile  = ",outputfile,"\n"))

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
  adataMatDT <- fread(inputfile, header=TRUE, sep="\t")
  adataMatDF <- as.data.frame(adataMatDT)
  rownames(adataMatDF) <- adataMatDF$V1; adataMatDF$V1 <- NULL
  # Source: https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
  # NOTE: The keyword assay is important to explicitely mention the counts
  adataSE <- SummarizedExperiment(assay=as.matrix(adataMatDF))

  # Annotate it with SingleR
  MouseRNAseqRefadataSE <- SingleR(test = adataSE, ref = MouseRNAseqRef, assay.type.test=1, labels = MouseRNAseqRef$label.main)
  ImmGenRefadataSE      <- SingleR(test = adataSE, ref = ImmGenRef, assay.type.test=1, labels = ImmGenRef$label.main)

  # Get output file names
  # outputImmGenAnnFile   <- "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/tregCNS/individual_samples/S509/data/05_markerGenes_S509_SingleR_ImmGenRef.txt"
  # outputMouseRnaAnnFile <- "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/tregCNS/individual_samples/S509/data/05_markerGenes_S509_SingleR_MouseRNAseqRef.txt"
	outputImmGenAnnFile   <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_SingleR_ImmGenRef.txt" ,sep='');
	outputMouseRnaAnnFile <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_SingleR_MouseRNAseqRef.txt" ,sep='');

  # Save output as text files
  write.table(data.frame("cellIDs" = rownames(ImmGenRefadataSE),ImmGenRefadataSE)          , file = outputImmGenAnnFile  , row.names = FALSE, col.names=T, sep = '\t', quote = F)
  write.table(data.frame("cellIDs" = rownames(MouseRNAseqRefadataSE),MouseRNAseqRefadataSE), file = outputMouseRnaAnnFile, row.names = FALSE, col.names=T, sep = '\t', quote = F)

}

##################### USER DEFINIED FUNCTIONS #########################
# Load Libraries
load_libraries <- function(){
	# Load libraries at start
	# if (!requireNamespace("BiocManager", quietly=TRUE))
  #     install.packages("BiocManager")
  # BiocManager::install("SingleR")
  suppressMessages(require(SingleR))
  suppressMessages(require(SummarizedExperiment))
  suppressMessages(require(data.table))
  suppressMessages(require(tools))
}

check_options <- function(){
	# Description of the script
	desc <- sprintf("
	----------------- SAMPLE USAGE ------------------
	- Rscript scripts/R_annotate_cells_using_singleR.R -if=/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/tregCNS/individual_samples/S509/data/04_markerGenes_S509_adata_T.matrix -of=/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/tregCNS/individual_samples/S509/data/05_markerGenes_S509.txt
	-------------------------------------------------
	CONTACT: 
		Gaurav Jain
		gaurav.jain@tum.de
	-------------------------------------------------\n
	")
	# create parser object
	parser <- ArgumentParser(description=cat(desc))

	# Add arguments 
	parser$add_argument("-if", "--inputfile"  , dest="inputfile"  , help="*Input mac2 peaks tab file"   , type="character", required=TRUE)
	parser$add_argument("-of", "--outputfile" , dest="outputfile" , help="*Output annotation file"      , type="character", required=TRUE)

	# Print the help message only if no arguments are supplied
	if(length(commandArgs(TRUE))==0){
		cat(desc)
		parser$print_help()
		quit()
	}

	# get command line options, if help option encountered print help and exit,
	# otherwise if options not found on command line then set defaults,
	args <- parser$parse_args()
	return(args)
}

## Call the main function in the end
main()

