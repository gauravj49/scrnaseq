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

# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()
%load_ext rpy2.ipython

# # Loading R libraries
%%R 
# library(SCnorm)
library(scran)
library(RColorBrewer)
# library(slingshot)
# library(monocle)
library(gam)
# library(clusterExperiment)
library(ggplot2)
library(plyr)
library(MAST)

# Reset random state
np.random.seed(2105)

#***************************************************
# Data processing steps
# 1) Reading and processing the input data
# 2) Quality control
# 3) Normalization + log transformation
# 4) Technical correction
# 5) Biological correction
# 6) Expression recovery (denoising)
#***************************************************
# Source: https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb
#***************************************************

# System variables and directories
projName        = "manec_tissues_all_samples_dca_first" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/tissue/01_preprocessing/{0}".format(projName); create_dir("{0}".format(output_dir))
cc_genes_file   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/Macosko_cell_cycle_genes.txt"
minGenesPerCell = 100
minCellsPergene = 2
bname           = projName
qcDir           = "{0}/qc".format(output_dir); create_dir(qcDir)
countsDir       = "{0}/counts".format(output_dir); create_dir(countsDir)

# 1 Reading the data
# https://github.com/theislab/scanpy/issues/267
tissueFilenames = [
                    'input/manec/pilot2/bulk1001_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/pilot2/bulk997_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/pilot2/bulk1018_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/pilot2/stomach1001_mouse_filtered_feature_bc_matrix.h5'
                  ]
adatas          = [sc.read_10x_h5(f) for f in tissueFilenames]

# Get the dictionary of tissue ids
tissueIdDict =  {
                    '0':'bulk1001', 
                    '1':'bulk997', 
                    '2':'bulk1018', 
                    '3':'stomach1001',
                    '4':'bulk1079'
                }

# Make the variable names unique for each annData object separately and run DCA on them
for i, adata in enumerate(adatas):
    tid = tissueIdDict[str(i)]
    adata.var_names_make_unique()
    # Perform basic prefiltering
    print("- {0} with shape {1}".format(tid, adata.to_df().shape))
    sc.pp.filter_cells(adata, min_genes=minGenesPerCell) 
    sc.pp.filter_genes(adata, min_cells=minCellsPergene)
    adata.var_names_make_unique()
    print(adata.to_df().shape)
    # Expression recovery (denoising) the data on raw counts
    sce.pp.dca(adata)
    # Save the denoised matrices
    adata.write("{0}/02_dca_{2}_{1}_matrix.h5ad".format(countsDir, projName, tid))

# Merge denoised datasets for different mices
adata = adatas[0].concatenate(adatas[1:])

# Make variable names unique
adata.var_names_make_unique()

# Add tissue id column for the batches
adata.obs['tissueID'] = adata.obs['batch'].map(tissueIdDict)

# Convert the sparse count matrices to dense represntation
adata.X = adata.X.toarray()

# AnnData object with n_obs × n_vars = 28603 × 12384 
#     obs: 'DCA_split', 'batch', 'n_counts', 'n_genes', 'size_factors', 'tissueID'
#     var: 'gene_ids', 'feature_types', 'genome', 'n_cells-0', 'n_cells-1', 'n_cells-2', 'n_cells-3', 'n_cells-4'


# Checking the total size of the data set
adata.shape # We have 28603 cells and 12384 genes in the dataset

# 2) Quality control 
# 2.1) Calculate QC covariates
adata.obs['n_counts'] = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes'] = (adata.X > 0).sum(1)             # This will be the same number in each cell as they are all imputed values for every gene

mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']

# 2.2) Plot QC metrics
# Sample quality plots
t1 = sc.pl.violin(adata, 'n_counts', groupby='tissueID', size=2, log=True, cut=0, show=False)
plt.savefig("{0}/02_dca_{1}_tissueID_nCounts_plot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
t2 = sc.pl.violin(adata, 'mt_frac', groupby='tissueID', show=False)
plt.savefig("{0}/02_dca_{1}_tissueID_mtFraction_plot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# 2.3) Data quality summary plots
p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/02_dca_{1}_genes_counts_scatterplot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

p2 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/02_dca_{1}_genes_counts_scatterplot_zoomedin.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# 2.3) Thresholding decision based on counts
p3 = sns.distplot(adata.obs['n_counts'], kde=False,bins=1000); #plt.show()
plt.savefig("{0}/02_dca_{1}_ncounts_histogramplot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
p4 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<2000], kde=False, bins=1000); #plt.show()
plt.savefig("{0}/02_dca_{1}_ncounts_histogramplot_lessthan_2000.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
p5 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>5000], kde=False, bins=1000); #plt.show()
plt.savefig("{0}/02_dca_{1}_ncounts_histogramplot_greaterthan_5000.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# 2.4) Thresholding decision based on genes
p6 = sns.distplot(adata.obs['n_genes'], kde=False, bins=1000); # plt.show()
plt.savefig("{0}/02_dca_{1}_genes_histogramplot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
p7 = sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, bins=500); # plt.show()
plt.savefig("{0}/02_dca_{1}_genes_histogramplot_lessthan_1000.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# 2.5) Filter cells according to identified QC thresholds:
origadata = adata.copy()
print('Total number of cells: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_counts = 700)
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

# 2.7) Save the filtered raw data as tab separated matrix 
# adata.to_df().to_csv("{0}/01_raw_T_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="CellId")
# adata.to_df().T.to_csv("{0}/01_raw_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol")

# 2.8) Calculations for the visualizations
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(adata, random_state = 2105)
sc.tl.umap(adata, random_state = 2105)

# 2.9) Plot visualizations
sc.pl.pca_scatter(adata, color='n_counts', show=False)
plt.savefig("{0}/02_dca_{1}_clustering_ncounts_PCA.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['tissueID'], show=False)
plt.savefig("{0}/02_dca_{1}_clustering_tissueID_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

########################
# 3) Normalization and log transformation
# Normalization using SCRAN
# This method requires a coarse clustering input to improve size factor esimation performance. 
# Thus, we use a simple preprocessing approach and cluster the data at a low resolution to get an 
# input for the size factor estimation. The basic preprocessing includes assuming all size factors 
# are equal (library size normalization to counts per million - CPM) and log-transforming the count data.

# Perform a clustering for scran normalization in clusters
adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)

# Preprocess variables for scran normalization
input_groups = adata_pp.obs['groups']
data_mat = adata.X.T

# Run scran in R
%%R -i data_mat -i input_groups -o size_factors
size_factors = computeSumFactors(data_mat, clusters=input_groups, min.mean=0.1)

# Delete adata_pp
del adata_pp

# Visualize the estimated size factors
adata.obs['size_factors'] = size_factors
sc.pl.scatter(adata, 'size_factors', 'n_counts', show=False)
plt.savefig("{0}/02_dca_{1}_scrna_sizefactors_vs_ncounts.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.scatter(adata, 'size_factors', 'n_genes' , show=False)
plt.savefig("{0}/02_dca_{1}_scrna_sizefactors_vs_ngenes.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sns.distplot(size_factors, bins=50, kde=False)
plt.savefig("{0}/02_dca_{1}_scrna_sizefactors_histogram.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Keep the count data in a counts layer
adata.layers["counts"] = adata.X.copy()

# Normalize adata 
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)

# Store the full data set in 'raw' as log-normalised data for statistical testing
adata.raw = adata

# 4) Technical correction
# 4.1) Batch Correction using ComBat
sc.pp.combat(adata, key='tissueID')

# 4.2) Highly Variable Genes
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig("{0}/03_dcaNorm_{1}_highly_variable_genes.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Highly variable gene information is stored automatically in the adata.var['highly_variable'] field. The dataset now contains:
# - A 'counts' layer with count data
# - log-normalized data in adata.raw
# - batch corrected data in adata.X
# - highly variable gene annotations in adata.var['highly_variable']
#   - The HVG labels will be used to subselect genes for clustering and trajectory analysis.

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
sc.pl.umap(adata, color=['S_score', 'G2M_score'], use_raw=False, show=False)
plt.savefig("{0}/03_dcaNorm{1}_S_G2M_Phase_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color='phase', use_raw=False, show=False)
plt.savefig("{0}/03_dcaNorm{1}_combined_Phase_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Save the normalized, log transformed, batch and cell cycle corrected data
adata.write("{0}/03_dcaNorm_{1}_matrix.h5ad".format(countsDir, projName))

# 7) Clustering
# 7.1) Perform clustering - using highly variable genes
sc.tl.louvain(adata, key_added='louvain', random_state=2105)
sc.tl.louvain(adata, resolution=1.0, key_added='louvain_r1'  , random_state=2105)
sc.tl.louvain(adata, resolution=0.9, key_added='louvain_r0.9', random_state=2105)
sc.tl.louvain(adata, resolution=0.8, key_added='louvain_r0.8', random_state=2105)
sc.tl.louvain(adata, resolution=0.7, key_added='louvain_r0.7', random_state=2105)
sc.tl.louvain(adata, resolution=0.6, key_added='louvain_r0.6', random_state=2105)
sc.tl.louvain(adata, resolution=0.5, key_added='louvain_r0.5', random_state=2105)
sc.tl.louvain(adata, resolution=0.4, key_added='louvain_r0.4', random_state=2105)
sc.tl.louvain(adata, resolution=0.3, key_added='louvain_r0.3', random_state=2105)
sc.tl.louvain(adata, resolution=0.2, key_added='louvain_r0.2', random_state=2105)
sc.tl.louvain(adata, resolution=0.1, key_added='louvain_r0.1', random_state=2105)

# Number of cells in each cluster
adata.obs['louvain'].value_counts()

for i in np.linspace(0.1,0.9,9):
  print(adata.obs['louvain_r{0:0.1f}'.format(i)].value_counts())
adata.obs['louvain_r1'].value_counts()

# 4.3) Visualizations
# Calculations for the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata, random_state = 2105)

# Plot visualizations
# Visualize the clustering and how this is reflected by different technical covariates
plt.figure(figsize=(100,8))
sc.pl.umap(adata, color=['louvain', 'louvain_r1', 'louvain_r0.9', 'louvain_r0.8', 'louvain_r0.7', 'louvain_r0.6', 'louvain_r0.5', 'louvain_r0.4', 'louvain_r0.3', 'louvain_r0.2', 'louvain_r0.1'], palette=sc.pl.palettes.default_64, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_dcaNorm_{1}_clustering_louvain_all_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# UMAPs
sc.pl.pca_scatter(adata, color='n_counts', size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_dcaNorm_{1}_clustering_ncounts_PCA.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['tissueID'], size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_dcaNorm_{1}_clustering_tissueID_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color='n_counts', size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_dcaNorm_{1}_clustering_ncounts_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['log_counts', 'mt_frac'], size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_dcaNorm_{1}_clustering_logCounts_mtFrac_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

