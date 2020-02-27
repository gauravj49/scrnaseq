# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)

# Loading the python libraries
import scanpy as sc
import scanpy.external as sce
from gjainPyLib import *
import pickle

import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
import anndata2ri

# USER DEFINED SETTINGS:
# Ignore R warning messages
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

# Reset random state
np.random.seed(2105)

# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()
%load_ext rpy2.ipython

# # Loading R libraries
%%R 
library(scran)
library(RColorBrewer)
library(gam)
library(ggplot2)
library(plyr)
library(MAST)

#***************************************************
# Data processing steps
# 1) Reading and processing the input data
# 2) Quality control
# 3) Normalization + log transformation
# 4) Technical correction
# 5) Biological correction
# 6) Expression recovery (denoising)
# 7) Feature selection
# 8) Dimensionality reduction
#***************************************************
# Source: https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb
#***************************************************

# System variables and directories
projName        = "manec_tissues_merged_except1079" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/tissue/01_mitTumorIdentification/{0}".format(projName); create_dir("{0}".format(output_dir))
cc_genes_file   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/Macosko_cell_cycle_genes.txt"
minGenesPerCell = 100
minCountPerCell = 300
minCellsPergene = 25
bname           = projName
qcDir           = "{0}/analysis".format(output_dir); create_dir(qcDir)
countsDir       = "{0}/data".format(output_dir); create_dir(countsDir)

# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

# 1 Reading the data
# https://github.com/theislab/scanpy/issues/267
tissueFilenames = [
                    'input/manec/tissue/bulk997_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/bulk1001_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/bulk1018_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/stomach1001_mouse_filtered_feature_bc_matrix.h5'
                  ]
adatas          = [sc.read_10x_h5(f) for f in tissueFilenames]

# Get the dictionary of tissue ids
tissueIdDict =  {
                    '0':'bulk1001', 
                    '1':'bulk997', 
                    '2':'bulk1018', 
                    '3':'stomach1001'
                }

# Make the variable names unique for each annData object separately and run DCA on them
for i, adata in enumerate(adatas):
    tid = tissueIdDict[str(i)]
    adata.var_names_make_unique()

    # Perform basic prefiltering
    print("- {0} with shape {1}".format(tid, adata.to_df().shape))
    sc.pp.filter_cells(adata, min_genes=1) 
    sc.pp.filter_genes(adata, min_cells=10)
    adata.var_names_make_unique()
    print(adata.to_df().shape)

