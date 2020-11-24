#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: findMotifs_in_rnasrqGenesList.py
- CONTACT: Gaurav Jain(gaurav.jain@dzne.edu)
***********************************************
"""
print (__doc__)

import argparse
import os.path
import sys
import textwrap
import string
import re
import numpy  as np
import matplotlib as mpl
#mp.use('Agg') # to use matplotlib without X11
import matplotlib.pyplot as plt
import math
import time
import commands
from gjainLIB import *
from numpy import loadtxt, dtype,float32
from termcolor import colored
from collections import *
from scipy import stats
import numpy.ma as ma
################ USER CONFIGURATION ###################
#######################################################

def main():
    # Get input options
    args         = check_options()
    input_file   = args.input_file
    output_file  = args.output_file
    l2fc         = args.l2fc
    padj         = args.padj
    pval         = args.pval
    baseMean     = args.base_mean
    
    # Create output dir if not exists
    create_dir(get_file_info(output_file)[0])

    commonGeneNames = 1
    if input_file.endswith("DE_RESULTS.txt"):
        commonGeneNames = 0
    elif input_file.endswith("DE_RESULTS_geneSymbols*.txt"):
        commonGeneNames = 1

    # Get the upregulated and downregulated files
    uprout_file = "{0}_UPRregulated_significant_genes.txt".format(get_file_info(output_file)[3])
    dwnout_file = "{0}_DWNregulated_significant_genes.txt".format(get_file_info(output_file)[3])
    output_file = "{0}_ALLregulated_significant_genes.txt".format(get_file_info(output_file)[3])

    # Get the output file handle
    fo = open(output_file, 'w')
    uo = open(uprout_file, 'w')
    do = open(dwnout_file, 'w')

    with open(input_file,'rU') as fg:
        # head Dm_WToffspring_APPfather_WTmother_over_Em_WToffspring_WTfather_APPmother_DE_RESULTS_common_gene_names.txt | cut -f1-9
        # GeneName	feature	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	Sample_1125DG_mm_rna_sr_magda_E_1_S16_L003_R1_001
        # Vma21	ENSMUSG00000073131	243.020533921707	-0.395115054037802	0.0934754621196951	-4.22693876101792	2.36892105207604e-05	0.0786610577022454	358.504139424987
        # Arhgef11	ENSMUSG00000041977	753.186740978746	0.321167415621671	0.0813581573133405	3.94757484962123	7.89468399031082e-05	0.0786610577022454	720.382435456327
        # Mtf2	ENSMUSG00000029267	569.133862045754	-0.383066672707145	0.0972216240875122	-3.94013858853391	8.1434550204728e-05	0.0786610577022454	743.157992549208
        # ....

        # Get the header and save it
        header = fg.readline().strip()
        
        # Initialize vars 
        skip_iv        = 0
        total_genes    = 0
        total_uprgenes = 0
        total_dwngenes = 0

        row_header=list()
        if commonGeneNames:
            # 0:GeneName	1:feature	2:baseMean	3:log2FoldChange	4:lfcSE	5:stat	6:pvalue	7:padj	8:Sample_9097_mm_rna_sr_Raaz_Y_1_S3_L001_R1_001
            column_header = header.split("\t")[8:]
        else:
            # 0:feature	1:baseMean	2:log2FoldChange	3:lfcSE	4:stat	5:pvalue	6:padj	7:Sample_9097_mm_rna_sr_Raaz_Y_1_S3_L001_R1_001
            column_header = header.split("\t")[9:]

        # Loop through rest of the file
        for line in useful_lines(fg):
            featureid        = str(line.split("\t")[0])
            if commonGeneNames:
                # 0:GeneName	1:feature	2:baseMean	3:log2FoldChange	4:lfcSE	5:stat	6:pvalue	7:padj	8:Sample_9097_mm_rna_sr_Raaz_Y_1_S3_L001_R1_001
                lbaseMean         = line.split("\t")[2]
                llog2FoldChange   = line.split("\t")[3]
                lpval             = line.split("\t")[6]
                lpadj             = line.split("\t")[7]
            else:
                # 0:feature	1:baseMean	2:log2FoldChange	3:lfcSE	4:stat	5:pvalue	6:padj	7:Sample_9097_mm_rna_sr_Raaz_Y_1_S3_L001_R1_001
                lbaseMean         = line.split("\t")[1]
                llog2FoldChange   = line.split("\t")[2]
                lpval             = line.split("\t")[5]
                lpadj             = line.split("\t")[6]

            # Skip genes with invalid log fold change i.e. "NA"
            if llog2FoldChange == 'NA' or lbaseMean == 'NA' or lpval == 'NA' or lpadj == 'NA':
                skip_iv += 1
                continue
            llog2FoldChange = float(llog2FoldChange)

            # Skip mirs with log2 fold change that are between (-0.5, 0.5)
            if -1*l2fc < llog2FoldChange < l2fc:
                skip_iv += 1
                continue

            # Filter for baseMean < 100 
            if float(lbaseMean) < baseMean:
                skip_iv += 1
                continue

            # If pvalue is set then skip the padj
            if pval != 1.00:
                if float(lpval) > pval:
                    skip_iv += 1
                    continue
            else:
                # Filter for padjusted > 0.05
                if float(lpadj) > padj:
                    skip_iv += 1
                    continue

            # Write the filtered genes to the output file
            fo.write("{0}\n".format(featureid))
            total_genes += 1

            if llog2FoldChange >= l2fc:
                uo.write("{0}\n".format(featureid))
                total_uprgenes += 1
            elif llog2FoldChange <= -1*l2fc:
                do.write("{0}\n".format(featureid))
                total_dwngenes += 1
            
        if total_genes > 0:
            # Print the filtered data summary
            print "\n- Conditions used:\n\t- log2FC    <> {0}\n\t- basemean  >= {1}\n\t- pvalue    <= {2}\n\t- padjusted <= {3}\n".format(l2fc, baseMean, pval, padj)
            print "\n- Total upregulated    features: {0}". format(total_uprgenes)
            print "- Total downregulated  features: {0}"  . format(total_dwngenes)
            print "- Total                features: {0}\n". format(total_genes)
        else:
            print "\n *********** NOTE **********"
            print "- Threshold conditions used:\n\t- log2FC    <> {0}\n\t- basemean  >= {1}\n\t- pvalue    <= {2}\n\t- padjusted <= {3}\n".format(l2fc, baseMean, pval, padj)
            print "\n- There is no data left after applying the filtering conditions. Your cut off values are too stringent. Please check threshold."
            print " ***************************"

    # Close the output file handles
    fo.close()
    uo.close()
    do.close()

    # Sorting the output file
    print "\n- Sorting the output file: {}".format(output_file)
    os.system("sort -u {0} -o {0}".format(output_file))
    print "\n- Sorting the output file: {}".format(uprout_file)
    os.system("sort -u {0} -o {0}".format(uprout_file))
    print "\n- Sorting the output file: {}".format(dwnout_file)
    os.system("sort -u {0} -o {0}".format(dwnout_file))

def useful_lines(f):
    ''' Filter out useless lines from the blast output file '''
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            continue
        if not line:
            continue
        yield line

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/findMotifs_in_rnaseqGenesList.py -if=output/christian/hypoxia1/02.1_christian_rnaseq_hypoxia1_ohneOutliers_rnaseq_bm25/results_DEseq2/HypIII4h_over_HypIIIK_DE_RESULTS_common_gene_names.txt -of=output/christian/hypoxia1/02.1_christian_rnaseq_hypoxia1_ohneOutliers_rnaseq_bm25/motifs/HypIII4h_over_HypIIIK_DE_selected_genesList.txt
        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gaurav.jain@dzne.de
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-if", metavar="--ifile" , help="*SmallRNA DEseq2 file", dest="input_file"  , type=str, required=True)
    parser.add_argument("-of", metavar="--ofile" , help="*Output file name"    , dest="output_file" , type=str, required=True)
    parser.add_argument("-lf", metavar="--l2fc"  , help=" Log2 Fold change (Default <> 0.5 )", dest="l2fc", type=float, default=0.5)
    parser.add_argument("-pj", metavar="--padj"  , help=" Adjusted pvalue  (Default <= 0.05)", dest="padj", type=float, default=0.05)
    parser.add_argument("-pv", metavar="--pval"  , help=" Normal   pvalue  (Default <= 1.00)", dest="pval", type=float, default=1.00)
    parser.add_argument("-bm", metavar="--bmean" , help=" Base Mean        (Default >= 100 )", dest="base_mean", type=int, default=100)

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

# main function
if __name__=="__main__":
      main()
