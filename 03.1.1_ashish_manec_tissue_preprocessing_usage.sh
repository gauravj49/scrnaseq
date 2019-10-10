# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)

import scanpy as sc
from gjainPyLib import *

# Merge 10x datasets for different mices
# https://github.com/theislab/scanpy/issues/267
# Source: https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb
tissueFilenames = [
                    'input/manec/pilot2/bulk1001_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/pilot2/bulk997_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/pilot2/bulk1018_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/pilot2/stomach1001_mouse_filtered_feature_bc_matrix.h5'
                  ]
adatas          = [sc.read_10x_h5(f) for f in tissueFilenames]
adata           = adatas[0].concatenate(adatas[1:])

# System variables and directories
projName        = "manec_tissues_merged_except1079" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/tissue/01_preprocessing/{0}".format(projName); create_dir("{0}".format(output_dir))
minGenesPerCell = 100
minCellsPergene = 2
bname = projName
qcDir = "{}/qc".format(projName); create_dir(qcDir)


# Add tissue id column for the batches
tissueIdDict =  {
                    '0':'bulk1001', 
                    '1':'bulk997', 
                    '2':'bulk1018', 
                    '3':'stomach1001'
                }
adata.obs['tissueID'] = adata.obs['batch'].map(tissueIdDict)

# Make variable names unique
adata.var_names_make_unique(); adata
# AnnData object with n_obs × n_vars = 9038 × 31053 
#     obs: 'batch', 'tissueID'
#     var: 'gene_ids', 'feature_types', 'genome'

# Checking the total size of the data set
adata.shape # We have 9038 cells and 31053 genes in the dataset

#***************************************************
# Data processing steps
# 1) Quality control
# 2) Normalization + log transformation
# 3) Technical correction
# 4) Biological correction
# 5) Feature selection
# 6) Dimensionality reduction
#***************************************************

# 1.1) Quality control - calculate QC covariates
adata.obs['n_counts'] = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes'] = (adata.X > 0).sum(1)

mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']

# 1.2) Quality control - plot QC metrics
# Sample quality plots
t1 = sc.pl.violin(adata, 'n_counts', groupby='tissueID', size=2, log=True, cut=0, show=False)
plt.savefig("{0}/{1}_tissueID_nCounts_plot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
t2 = sc.pl.violin(adata, 'mt_frac', groupby='tissueID', show=False)
plt.savefig("{0}/{1}_tissueID_mtFraction_plot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')