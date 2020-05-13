#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: deseq2_comparisons_venn.py
- CONTACT: Gaurav Jain(gaurav.jain@dzne.edu)
***********************************************
"""

################ USER CONFIGURATION ###################

#######################################################

def main():
    # Print initial stats
    nFiles = len(input_files_list) 
    print("- Number of files: {0}" .format(nFiles))

    # Define the data frame of features with file names as column names
    upfdf = pd.DataFrame()

    # Loop through the files and get the common features
    print("- Looping through the input files and getting the common features ...")
    for i, input_file in enumerate(input_files_list):
    
        print("\t- Processing: {0} ...".format(get_file_info(input_file)[4]))
        # 1) Get the differentially expressed features as a list
        upff = get_significant_features(input_file, lfc, padj, pval, score)

        # Get the column name
        if fileLabels:
            columnName = fileLabels[i]
        else:
            columnName = get_file_info(input_file)[1]

        # 2) Convert the list to pandas dataframe and concatenate the dataframe with "df" above
        upfdf = pd.concat([upfdf, pd.DataFrame(upff, columns=[columnName])], axis=1)


    # Save and draw output files for significant marker genes
    print("\n- Saving and drawing output files for significant marker genes") 
    save_draw_venn(upfdf, "{0}_significant_marker_genes.txt".format(get_file_info(output_file)[3])  , venn_plots_dir)

################ USER DEFINED FUNCTIONS ###################
def save_draw_venn(df, output_file, venn_plots_dir):
    ''' Save the output file and draw the venns'''

    # Save the dataframe to output file
    print("\n{0}".format('#' * 50))
    print("- Saving dataframe to output file: {0}".format(output_file))
    df.to_csv(output_file, sep='\t', index=False)

    # Draw the venn diagrams
    os.system("python {0}/scripts/venn_plots.py -if={1} -of={2}".format(venn_plots_dir, output_file, "{0}.png".format(get_file_info(output_file)[3])))
    

def get_significant_features(input_file, lfc_param, padj_param, pval_param, score_param):
    ''' Parse the input file and get the log2foldchange for each gene '''

    # Read the file in a dataframe
    # names   scores  logfoldchanges  pvals   pvals_adj
    # Ptma    8.1     1.1     4.7e-15 5.8e-12
    # Npm1    7.3     1.1     1.1e-12 6.2e-10
    # Ppia    7.1     1       4e-12   1.5e-09
    # Eef1a1  6.1     0.73    2e-09   3.1e-07
    # Hspd1   5.9     1.1     6.5e-09 8.4e-07
    rankedGenesDF = pd.read_csv(input_file, sep='\t', quoting=3, index_col=0)
    if score_param:
        significantMarkerGenesDF = rankedGenesDF[rankedGenesDF['scores'] >= score_param]
    else:
        # condition for logfc, pvalue or padj
        significantMarkerGenesDF = rankedGenesDF[((rankedGenesDF['logfoldchanges'] <= -1*lfc_param)|((rankedGenesDF['logfoldchanges'] >= lfc_param))) & (rankedGenesDF['pvals_adj'] <= padj_param)]

    # Get sorted unique list
    significantMGlist = sorted(get_unique_list(significantMarkerGenesDF.index.values.tolist()))
    print("\t\t- Number of significant marker genes: {0}".format(len(significantMGlist)))

    return list(significantMGlist)

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/marker_list_comparisons_venn.py -if "output/manec/individual_samples/bulk997/data/rankedGenes/louvain_r05/03_bulk997_rank_genes_louvain_r05_0.txt" "output/manec/individual_samples/bulk1001/data/rankedGenes/louvain_r1/03_bulk1001_rank_genes_louvain_r1_0.txt" "output/manec/individual_samples/bulk1018/data/rankedGenes/louvain_r1/03_bulk1018_rank_genes_louvain_r1_0.txt" -of=test.txt -sc=5
        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gaurav.jain@dzne.de
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-if", metavar="--inflst", help="*List of files to compapre" , dest="input_files_list", type=str, nargs='+', required=True)
    parser.add_argument("-of", metavar="--ofile" , help="*Output file name"          , dest="output_file" , type=str, required=True)
    parser.add_argument("-fl", metavar="--filab" , help=" File labels to use"        , dest="fileLabels"  , type=str)
    parser.add_argument("-lf", metavar="--lfc"   , help=" Log2 Fold change (Default <> 0.5 )", dest="lfc" , type=float, default=0.5)
    parser.add_argument("-pj", metavar="--padj"  , help=" Adjusted pvalue  (Default <= 0.05)", dest="padj", type=float, default=0.05)
    parser.add_argument("-pv", metavar="--pval"  , help=" Normal   pvalue  (Default <= 1.00)", dest="pval", type=float, default=1.00)
    parser.add_argument("-sc", metavar="--score" , help=" z-score underlying the computation of a p-value for each gene for each group\n\t- If scores is set then skip logfoldchanges and padjusted thresholds"  , dest="score",type=float)
    parser.add_argument("-vd", metavar="--vndir" , help=" Venn plots dir: (Default /home/rad/packages/venn_plots )", dest="venn_plots_dir", type=str, default="/home/rad/packages/venn_plots")

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().output_file:
        create_dir(get_file_info(parser.parse_args().output_file)[0])
        logfile = get_file_info(parser.parse_args().output_file)[3] + ".log"
    else:
        cwd = os.getcwd()
        logfile = "{0}/{1}.log".format(cwd,get_file_info(sys.argv[0])[1])
    logf = open(logfile, 'w')
    sys.stdout = Log(logf, sys.stdout)

    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()
    print_initial_arguments(parser)
    return args

if __name__=="__main__":
    print (__doc__)

    # Built in modules
    import argparse
    import os.path
    import sys

    # 3rd party modules
    import textwrap
    import re
    import numpy as np
    import scipy as sp
    import pandas as pd
    import matplotlib as mp
    #mp.use('Agg') # to use matplotlib without X11
    import matplotlib.pyplot as plt
    import subprocess
    import binascii as bi
    import scipy.stats as stats
    from collections import *
    from numpy import nanmean

    # for looping files in a dir
    import glob

    # user defined modules
    from gjainPyLib import *      # import all the functions from the Gaurav`s python library

    ### for color scale
    from  matplotlib import colors
    from itertools import cycle, islice # barplot colors

    ################ USER CONFIGURATION ###################
    np.set_printoptions(precision=6)
    #######################################################

    # Get input options
    args = check_options()

    # Store the variables
    input_files_list = [f for f in args.input_files_list]
    output_file      = args.output_file
    lfc              = args.lfc
    padj             = args.padj
    pval             = args.pval
    score            = args.score
    fileLabels       = args.fileLabels
    venn_plots_dir   = args.venn_plots_dir
    
    # Get output file name
    ofile = open(output_file, 'w')

    # Get the filelables in the list
    if fileLabels:
        fileLabels = fileLabels.split(",")


    main()

