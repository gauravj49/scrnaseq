# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)

# Loading the python libraries
import scanpy as sc
from gjainPyLib import *
import pickle

import rpy2.rinterface_lib.callbacks
import logging
from rpy2.robjects import pandas2ri
import anndata2ri

# USER DEFINED SETTINGS:
# Ignore R warning messages
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()
%load_ext rpy2.ipython

# Loading R libraries
%%R 
library(SCnorm)

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

# System variables and directories
projName        = "manec_tissues_merged_except1079" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/tissue/01_preprocessing/{0}".format(projName); create_dir("{0}".format(output_dir))
cc_genes_file   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/Macosko_cell_cycle_genes.txt"
minGenesPerCell = 100
minCellsPergene = 2
bname           = projName
qcDir           = "{0}/qc".format(output_dir); create_dir(qcDir)
countsDir       = "{0}/counts".format(output_dir); create_dir(countsDir)

# 1 Reading the data
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

# Make variable names unique
adata.var_names_make_unique()

# Add tissue id column for the batches
tissueIdDict =  {
                    '0':'bulk1001', 
                    '1':'bulk997', 
                    '2':'bulk1018', 
                    '3':'stomach1001'
                }
adata.obs['tissueID'] = adata.obs['batch'].map(tissueIdDict)

# Convert the sparse count matrices to dense represntation
adata.X = adata.X.toarray()

# AnnData object with n_obs × n_vars = 9038 × 31053 
#     obs: 'batch', 'tissueID'
#     var: 'gene_ids', 'feature_types', 'genome'

# Checking the total size of the data set
adata.shape # We have 9038 cells and 31053 genes in the dataset

# 2) Quality control 
# 2.1) Calculate QC covariates
adata.obs['n_counts'] = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes'] = (adata.X > 0).sum(1)

mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']

# 2.2) Plot QC metrics
# Sample quality plots
t1 = sc.pl.violin(adata, 'n_counts', groupby='tissueID', size=2, log=True, cut=0, show=False)
plt.savefig("{0}/{1}_tissueID_nCounts_plot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
t2 = sc.pl.violin(adata, 'mt_frac', groupby='tissueID', show=False)
plt.savefig("{0}/{1}_tissueID_mtFraction_plot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# 2.3) Data quality summary plots
p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/{1}_genes_counts_scatterplot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

p2 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/{1}_genes_counts_scatterplot_zoomedin.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# NOTE: 
# It can be seen in the main cloud of data points, that cells with lower counts and genes tend to have a higher fraction of mitochondrial counts. 
# These cells are likely under stress or are dying. 
# When apoptotic cells are sequenced, there is less mRNA to be captured in the nucleus, and 
# therefore fewer counts overall, and thus a higher fraction of counts fall upon mitochondrial RNA. 
# If cells with high mitochondrial activity were found at higher counts/genes per cell, 
# this would indicate biologically relevant mitochondrial activity.

# 2.3) Thresholding decision based on counts
p3 = sns.distplot(adata.obs['n_counts'], kde=False); #plt.show()
plt.savefig("{0}/{1}_ncounts_histogramplot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
p4 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<2000], kde=False, bins=1000); #plt.show()
plt.savefig("{0}/{1}_ncounts_histogramplot_lessthan_2000.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
p5 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>5000], kde=False, bins=1000); #plt.show()
plt.savefig("{0}/{1}_ncounts_histogramplot_greaterthan_5000.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# 2.4) Thresholding decision based on genes
p6 = sns.distplot(adata.obs['n_genes'], kde=False, bins=1000); # plt.show()
plt.savefig("{0}/{1}_genes_histogramplot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
p7 = sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, bins=500); # plt.show()
plt.savefig("{0}/{1}_genes_histogramplot_lessthan_1000.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')


# 2.5) Filter cells according to identified QC thresholds:
origadata = adata.copy()
print('Total number of cells: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_counts = 300)
print('Number of cells after min count filter: {:d}'.format(adata.n_obs))

# sc.pp.filter_cells(adata, max_counts = 40000)
# print('Number of cells after max count filter: {:d}'.format(adata.n_obs))

adata = adata[adata.obs['mt_frac'] < 0.25]
print('Number of cells after MT filter: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_genes = 300)
print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

# Total number of cells: 9038
# Number of cells after min count filter: 9038
# Number of cells after MT filter: 7992
# Trying to set attribute `.obs` of view, making a copy.
# Number of cells after gene filter: 7261

# 2.6) Filter genes according to identified QC thresholds:
# Min 5 cells - filters out 0 count genes
print('Total number of genes: {:d}'.format(adata.n_vars))
sc.pp.filter_genes(adata, min_cells=5)
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))

# Total number of genes: 31053
# Number of genes after cell filter: 16263

# Save the filtered raw data as tab separated matrix 
adata.to_df().to_csv("{0}/01_raw_T_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="CellId")
adata.to_df().T.to_csv("{0}/01_raw_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol")

# 3) Normalization + log transformation

# Get the raw matrix to pass it along to Scnorm
rawMat       = adata.to_df().T
input_groups = adata.obs['tissueID'].tolist()

with open('conditions.txt', 'w') as f: f.write("\n".join(str(item) for item in input_groups))

# Run Scnorm in R
%%R -i rawMat -i input_groups -o normRawData
groups <- (as.vector(unlist(input_groups)))
normRawMat  <- SCnorm(Data=rawMat, Conditions=groups, K=2)
normRawData <- results(normRawMat, type = "NormalizedData")

# Keep the count data in a counts layer
adata.layers["counts"] = adata.X.copy()

# Normalize and log adata 
adata.X = normRawData
sc.pp.log1p(adata)

# 4) Technical correction
# 4.1) Batch Correction using ComBat
sc.pp.combat(adata, key='tissueID')

# 4.2) Highly Variable Genes
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))

# Highly variable gene information is stored automatically in the adata.var['highly_variable'] field. The dataset now contains:
# a 'counts' layer with count data
# log-normalized data in adata.raw
# batch corrected data in adata.X
# highly variable gene annotations in adata.var['highly_variable']
# The HVG labels will be used to subselect genes for clustering and trajectory analysis.


# 5) Biological correction
# 5.1) Cell cycle scoring
#Score cell cycle and visualize the effect:
cc_genes         = pd.read_table(cc_genes_file, delimiter='\t')
s_genes          = cc_genes['S'].dropna()
g2m_genes        = cc_genes['G2.M'].dropna()
s_genes_mm       = [gene.lower().capitalize() for gene in s_genes]
g2m_genes_mm     = [gene.lower().capitalize() for gene in g2m_genes]
s_genes_mm_ens   = adata.var_names[np.in1d(adata.var_names, s_genes_mm)]
g2m_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, g2m_genes_mm)]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)
sc.pl.umap(adata, color=['S_score', 'G2M_score'], use_raw=False)
plt.savefig("{0}/{1}_S_G2M_Phase_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color='phase', use_raw=False)
plt.savefig("{0}/{1}_combined_Phase_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Save the normalized, log transformed, batch and cell cycle corrected data
adata.to_df().to_csv("{0}/02_normalizedRaw_T_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="CellId")
adata.to_df().T.to_csv("{0}/02_normalizedRaw_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol")


# 6) Expression recovery (denoising)
# Denoise the data on raw counts
dca(adata)

# Save the denoised data
adata.to_df().to_csv("{0}/02_denoised_DCA_T_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="CellId")
adata.to_df().T.to_csv("{0}/02_denoised_DCA_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol")

# 7) Feature selection
# 8) Dimensionality reduction