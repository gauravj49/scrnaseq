
##################################################################
ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)

# Loading the python libraries
import scanpy as sc
import scanpy.external as sce
from gjainPyLib import *
import pickle
import logging
import scanorama
import trvae
import scvelo
from textwrap import wrap

# For X11 display
import matplotlib
# matplotlib.use('TkAgg')
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D      # for 3D projection
from matplotlib.colors import ListedColormap # for sc.palette to colormap
from itertools import combinations           # pairwise combinations

# Reset random state
np.random.seed(2105)

# For using R inside python
import rpy2.rinterface_lib.callbacks
from rpy2.robjects import pandas2ri
import anndata2ri

# Ignore R warning messages
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

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
# 4) Biological correction
# 5) Expression recovery (denoising)
# 6) Technical correction
# 7) Feature selection
# 8) Dimensionality reduction
#***************************************************
# Source: https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb
#***************************************************

# System variables and directories
projName        = "tregCNS"
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/{0}/allSamples".format(projName); create_dir("{0}".format(output_dir))
ccGenes_macosko = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/macosko_cell_cycle_genes_mmu.txt"
ccGenes_regev   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/regev_lab_cell_cycle_genes_mmu.txt"
minGenesPerCell = 5
minCountPerCell = 50
maxCountPerCell = 100000 
minCellsPergene = 2
mtGenesFilter   = 0.25
rbGenesFilter   = 0.30
bname           = projName
plotsDir        = "{0}/plots".format(output_dir); create_dir(plotsDir)
dataDir         = "{0}/data".format(output_dir); create_dir(dataDir)

# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

# Get the logging parameters
sc.logging.print_memory_usage()
sc.logging.print_version_and_date()
sc.logging.print_versions()
# Memory usage: current 1.07 GB, difference +1.07 GB
# Running Scanpy 1.4.6, on 2020-07-16 23:51.
# scanpy==1.4.6 anndata==0.7.3 umap==0.4.4 numpy==1.19.0rc2 scipy==1.4.1 pandas==0.25.3 scikit-learn==0.23.1 statsmodels==0.11.1 python-igraph==0.8.2 louvain==0.7.0

# 1) Reading and performing QC on individual datasets
# 1.1) Reading the data in the anndata object individually
# adata = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S509')
# S503/  S504/  S505/  S508/  S509/  S511/  S512/  S514/  S515/  S516/  S517/  S518/  S519/
tissueFilenames = [
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S503',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S504',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S505',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S508',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S509',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S511',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S512',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S514',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S515',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S516',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S517',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S518',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S519'
                  ]
adatas          = [sc.read_10x_mtx(f) for f in tissueFilenames]
adatas
# Out[8]:
# [AnnData object with n_obs × n_vars =   1925 × 55471
#  AnnData object with n_obs × n_vars =   3363 × 55471
#  AnnData object with n_obs × n_vars =   4516 × 55471
#  AnnData object with n_obs × n_vars =    137 × 55471
#  AnnData object with n_obs × n_vars =   1432 × 55471
#  AnnData object with n_obs × n_vars =    983 × 55471
#  AnnData object with n_obs × n_vars =    974 × 55471
#  AnnData object with n_obs × n_vars =   2944 × 55471
#  AnnData object with n_obs × n_vars = 147456 × 55471
#  AnnData object with n_obs × n_vars =   3789 × 55471
#  AnnData object with n_obs × n_vars =   2614 × 55471
#  AnnData object with n_obs × n_vars =   2187 × 55471
#  AnnData object with n_obs × n_vars =   3768 × 55471


# 1.2) Get the dictionary of tissue ids
tissueIdDict =  {
                    '0' :'S503', 
                    '1' :'S504', 
                    '2' :'S505', 
                    '3' :'S508', 
                    '4' :'S509', 
                    '5' :'S511', 
                    '6' :'S512', 
                    '7' :'S514', 
                    '8' :'S515', 
                    '9' :'S516', 
                    '10':'S517', 
                    '11':'S518', 
                    '12':'S519', 
                }
# adata3 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S503')
# adata4 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S504')
# adata5 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S505')
# adata8 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S508')
# adata9 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S509')
# adata11 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S511')
# adata12 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S512')
# adata14 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S514')
# adata15 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S515')
# adata16 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S516')
# adata17 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S517')
# adata18 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S518')
# adata19 = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S519')

# 1.3) Merge 10x datasets for different mices
# https://github.com/theislab/scanpy/issues/267
adata = adatas[0].concatenate(adatas[1:])

# 1.4) Make variable names unique
adata.var_names_make_unique()
# Convert the sparse count matrices to dense represntation
adata.X = adata.X.toarray()

# 1.5) Add tissue id column for the batches
adata.obs['tissueID'] = adata.obs['batch'].map(tissueIdDict)
# In [16]: adata.obs
# Out[16]:
#                     batch tissueID
# AAACAAACAGCTATGA-0      0     S503
# AAACAAACCTACGAGC-0      0     S503

# 1.6) Checking the total size of the data set
adata.shape # We have 5887 cells and 6068 genes in the dataset


# 1.2.1) Calculate QC covariates
print("- Shape {0}".format(adata.to_df().shape))
sc.pp.calculate_qc_metrics(adata, inplace=True) # we now have many additional data types in the obs slot:
adata.obs['n_counts']   = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes']    = (adata.X > 0).sum(1)
adata
# obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'n_counts', 'log_counts', 'n_genes'
# var: 'gene_ids', 'feature_types', 'genome', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'

# 1.2.2) Calculate mitochondrial/ribosomal genes fraction (percentage)
# For each cell compute fraction of counts in mito/ribo genes vs. all genes
mt_gene_mask         = [gene.startswith('mt-') for gene in adata.var_names]
adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']
rb_gene_mask         = [gene.startswith(("Rps","Rpl")) for gene in adata.var_names]
adata.obs['rb_frac'] = adata.X[:, rb_gene_mask].sum(1)/adata.obs['n_counts']

# 1.2.3) Plot QC metrics
fig = plt.figure(figsize=(15,20))
n   = 
# Sample quality plots
ax = fig.add_subplot(n, 2, 1); t1  = sc.pl.violin(adata, ['n_genes_by_counts', 'n_counts'], jitter=0.4, size=2, log=True, cut=0, ax = ax, show=False)
ax = fig.add_subplot(n, 2, 2); t2  = sc.pl.violin(adata, ['mt_frac','rb_frac'], jitter=0.4, size=2, log=False, cut=0, ax = ax, show=False)
# 1.2.4) Thresholdingecision based on counts
ax = fig.add_subplot(n, 2, 3); p3  = sns.distplot(adata.obs['n_counts'], kde=False, ax = ax, bins=50); #plt.show()
ax = fig.add_subplot(n, 2, 4); p4  = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<2000], kde=False, ax = ax, bins=50); #plt.show()
# 1.2.5) Thresholding decision based on genes
ax = fig.add_subplot(n, 2, 5); p6  = sns.distplot(adata.obs['n_genes'], kde=False, ax = ax, bins=50); # plt.show()
ax = fig.add_subplot(n, 2, 6); p7  = sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, ax = ax, bins=50); # plt.show()
# 1.2.6) mt and rb fraction plots
ax = fig.add_subplot(n, 2, 7); p8  = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', show=False)
ax = fig.add_subplot(n, 2, 8); p9  = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', show=False)
ax = fig.add_subplot(n, 2, 9); p10 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='rb_frac', show=False)
ax = fig.add_subplot(n, 2, 10);p11 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='rb_frac', show=False)
plt.tight_layout()
plt.savefig("{0}/01_raw_{1}_QC_matrices.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 1.2.6) Data quality summary plots
p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/01_raw_{1}_genes_counts_mtfrac_scatterplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p2 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/01_raw_{1}_genes_counts_mtfrac_scatterplot_zoomedin.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='rb_frac', show=False)
plt.savefig("{0}/01_raw_{1}_genes_counts_rbfrac_scatterplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p2 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='rb_frac', show=False)
plt.savefig("{0}/01_raw_{1}_genes_counts_rbfrac_scatterplot_zoomedin.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 1.2.7) Filter cells according to identified QC thresholds:
print('Total number of cells: {:d}'.format(adata.n_obs))
sc.pp.filter_cells(adata, min_counts = minCountPerCell)
print('Number of cells after min count filter: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, max_counts = maxCountPerCell)
print('Number of cells after max count filter: {:d}'.format(adata.n_obs))

adata = adata[adata.obs['mt_frac'] < mtGenesFilter]
print('Number of cells after MT filter  : {:d}'.format(adata.n_obs))
adata = adata[adata.obs['rb_frac'] < rbGenesFilter]
print('Number of cells after Ribo filter: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_genes = minGenesPerCell)
print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

# 1.2.8) Filter genes according to identified QC thresholds:
print('Total number of genes: {:d}'.format(adata.n_vars))
sc.pp.filter_genes(adata, min_cells=minCellsPergene)
print('Number of genes after minCellsPergene filter: {:d}'.format(adata.n_vars))

# Total number of cells                 : 176088
# Number of cells after min count filter: 58834
# Number of cells after max count filter: 58834
# Number of cells after MT filter       : 58782
# Number of cells after Ribo filter     : 57975
# Total number of genes: 55471
# Number of genes after minCellsPergene filter: 27119

# 1.2.9) Compute variable genes
# We first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.
sc.pp.highly_variable_genes(adata, flavor='cell_ranger')
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))

# 1.2.10) Calculations for the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(adata, random_state = 2105)
sc.tl.umap(adata, random_state = 2105, n_components=3)

# 1.2.11) Plot visualizations
sc.pl.pca_scatter(adata, color='n_counts',show=False)
plt.savefig("{0}/01_raw_{1}_clustering_ncounts_PCA.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
# UMAPS
fig = plt.figure(figsize=(20,12))
# 2D projection
ax = fig.add_subplot(2, 3, 1);                  sc.pl.umap(adata   ,                  ax=ax, color='log_counts'   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format('log_counts'))
ax = fig.add_subplot(2, 3, 2);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="mt_frac", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="mt_frac UMAP")
ax = fig.add_subplot(2, 3, 3);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="rb_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="rb_frac UMAP")
# 3D projection
ax = fig.add_subplot(2, 3, 4, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color='log_counts', palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format('log_counts'))
ax = fig.add_subplot(2, 3, 5, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="mt_frac UMAP")
ax = fig.add_subplot(2, 3, 6, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color="rb_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="rb_frac UMAP")
plt.tight_layout()
plt.savefig("{0}/01_raw_{1}_logCounts_mt_rb_frac_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# # 1.2.12) Expression recovery (denoising) the data on raw counts
# sce.pp.dca(adata)

# 1.3) Checking the total size of the data set
adata.shape # We have 5887 cells and 11029 genes in the dataset

# 1.4) Save the filtered raw adata into a file
# Write the adata object to file
adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the filtered raw adata object
# adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata  = sc.read_h5ad(adatafile)

########################
# Normalization using CPM
########################
# Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells
cpmadata = adata.copy()
sc.pp.normalize_total(cpmadata, target_sum=1e4)
sc.pp.log1p(cpmadata) # Logarithmize the data.
cpmadata.raw = cpmadata

# 3) # Normalization using SCRAN
scranadata = adata.copy()
# Perform a clustering for scran normalization in clusters
adata_pp = scranadata.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)
# Preprocess variables for scran normalization
input_groups = adata_pp.obs['groups']
data_mat = scranadata.X.T
# Run scran in R
%%R -i data_mat -i input_groups -o size_factors
size_factors = computeSumFactors(data_mat, clusters=input_groups, min.mean=0.1)
# Delete adata_pp
del adata_pp
# Visualize the estimated size factors
scranadata.obs['size_factors'] = size_factors
fig = plt.figure(figsize=(16,6))
fig.suptitle('Estimated size factors')
ax = fig.add_subplot(1, 2, 1)
sc.pl.scatter(scranadata, 'size_factors', 'n_counts', ax=ax, show=False)
ax = fig.add_subplot(1, 2, 2)
sc.pl.scatter(scranadata, 'size_factors', 'n_genes', ax=ax, show=False)
plt.tight_layout()
plt.savefig("{0}/02_norm_{1}_scran_sizefactors_plots.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
# Keep the count data in a counts layer
scranadata.layers["counts"] = scranadata.X.copy()
# Normalize adata 
scranadata.X /= scranadata.obs['size_factors'].values[:,None]
sc.pp.log1p(scranadata)
# Store the full data set in 'raw' as log-normalised data for statistical testing
scranadata.raw = scranadata

# 4) Biological correction
# 4.1) Read cell cycle genes
cc_genes         = pd.read_table(ccGenes_macosko, delimiter='\t')
s_genes          = cc_genes['S'].dropna()
g2m_genes        = cc_genes['G2.M'].dropna()
# For mouse only
s_genes_mm       = [gene.lower().capitalize() for gene in s_genes]
g2m_genes_mm     = [gene.lower().capitalize() for gene in g2m_genes]

# CPM 
s_genes_mm_ens   = cpmadata.var_names[np.in1d(cpmadata.var_names, s_genes_mm)]
g2m_genes_mm_ens = cpmadata.var_names[np.in1d(cpmadata.var_names, g2m_genes_mm)]
# SCRAN
s_genes_mm_ens   = scranadata.var_names[np.in1d(scranadata.var_names, s_genes_mm)]
g2m_genes_mm_ens = scranadata.var_names[np.in1d(scranadata.var_names, g2m_genes_mm)]

# 4.2) Score cell cycle genes
sc.tl.score_genes_cell_cycle(cpmadata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)
sc.tl.score_genes_cell_cycle(scranadata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)

# Compute variable genes
# We first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.
sc.pp.highly_variable_genes(scranadata, flavor='cell_ranger')
print('\n','Number of highly variable genes: {:d}'.format(np.sum(scranadata.var['highly_variable'])))

# Calculations for the visualizations
sc.pp.pca(scranadata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(scranadata)
sc.tl.umap(scranadata, random_state = 2105, n_components=3)

# 4.3) Visualize the effects of cell cycle
fig = plt.figure(figsize=(16,12))
fig.suptitle('Effects of Cell Cycle')
ax = fig.add_subplot(2, 2, 1)
sc.pl.umap(scranadata, color=['S_score']  , ax=ax, use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
ax = fig.add_subplot(2, 2, 2)
sc.pl.umap(scranadata, color=['G2M_score'], ax=ax, use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
ax = fig.add_subplot(2, 2, 3)
sc.pl.umap(scranadata, color='phase', ax=ax, use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
ax = fig.add_subplot(2, 2, 4, projection='3d')
sc.pl.umap(scranadata, color='phase', ax=ax, use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.tight_layout()
plt.savefig("{0}/02_norm_{1}_scran_cell_cycle_plots.png".format(plotsDir, bname) , bbox_inches='tight', dpi=75); plt.close('all')


# 4) Technical correction
adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata  = sc.read_h5ad(adatafile)
# 4.1) Batch Correction using Combat
# CPM
cpmcombatadata = cpmadata.copy()
sc.pp.combat(cpmcombatadata, key='tissueID')
sc.pp.highly_variable_genes(cpmcombatadata, flavor='cell_ranger', n_top_genes=4000)
sc.pp.pca(cpmcombatadata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(cpmcombatadata, random_state = 2105)
sc.tl.umap(cpmcombatadata, random_state = 2105, n_components=3)
#SCRAN
scrancombatadata = scranadata.copy()
sc.pp.combat(scrancombatadata, key='tissueID')
sc.pp.highly_variable_genes(scrancombatadata, flavor='cell_ranger', n_top_genes=4000)
sc.pp.pca(scrancombatadata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(scrancombatadata, random_state = 2105)
sc.tl.umap(scrancombatadata, random_state = 2105, n_components=3)

# 4.2.1) Batch Correction using Scanorama for cpm normalized data
cpmadata2 = sc.AnnData(X=cpmadata.X, var=cpmadata.var, obs = cpmadata.obs)
#variable genes for the full dataset
sc.pp.highly_variable_genes(cpmadata2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'batch')
var_genes_batch = cpmadata2.var.highly_variable_nbatches > 0
var_select = cpmadata2.var.highly_variable_nbatches > 1
var_genes = var_select.index[var_select]
# Split per batch into new objects.
batches = ['0','1','2','3','4','5','6','7','8','9','10','11','12']
cpmalldata = {}
for batch in batches:
    cpmalldata[batch] = cpmadata2[cpmadata2.obs['batch'] == batch,]
# Subset the individual dataset to the same variable genes as in MNN-correct.
cpmalldata2 = dict()
for ds in cpmalldata.keys():
    print(ds)
    cpmalldata2[ds] = cpmalldata[ds][:,var_genes]
# Convert to list of AnnData objects
comnormadatas = list(cpmalldata2.values())
# Run scanorama.integrate
cpmscanorama  = scanorama.integrate_scanpy(comnormadatas, dimred = 50,)
# Make into one matrix.
cpmall_s = np.concatenate(cpmscanorama)
print(cpmall_s.shape)
# Add to the AnnData object
cpmscanoramaadata = adata.copy()
cpmscanoramaadata.obsm["SC"] = cpmall_s
# Calculations for the visualizations
sc.pp.highly_variable_genes(cpmscanoramaadata, flavor='cell_ranger', n_top_genes=4000)
sc.pp.pca(cpmscanoramaadata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(cpmscanoramaadata, random_state = 2105, use_rep = "SC")
sc.tl.umap(cpmscanoramaadata, random_state = 2105, n_components=3)

# 4.2.2) Batch Correction using Scanorama for scrna normalized data
scranadata2 = sc.AnnData(X=scranadata.X, var=scranadata.var, obs = scranadata.obs)
#variable genes for the full dataset
sc.pp.highly_variable_genes(scranadata2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'batch')
var_genes_batch = scranadata2.var.highly_variable_nbatches > 0
var_select = scranadata2.var.highly_variable_nbatches > 1
var_genes = var_select.index[var_select]
# Split per batch into new objects.
# batches = ['0','1','2','3','4','5','6','7','8','9','10','11','12']
batches = list(tissueIdDict.keys())
scranalldata = {}
for batch in batches:
    scranalldata[batch] = scranadata2[scranadata2.obs['batch'] == batch,]
# Subset the individual dataset to the same variable genes as in MNN-correct.
scranalldata2 = dict()
for ds in scranalldata.keys():
    print(ds)
    scranalldata2[ds] = scranalldata[ds][:,var_genes]
# Convert to list of AnnData objects
comnormadatas = list(scranalldata2.values())
# Run scanorama.integrate
scranscanorama  = scanorama.integrate_scanpy(comnormadatas, dimred = 50,)
# Make into one matrix.
scranall_s = np.concatenate(scranscanorama)
print(scranall_s.shape)
# Add to the AnnData object
scranscanoramaadata = adata.copy()
scranscanoramaadata.obsm["SC"] = scranall_s
# Calculations for the visualizations
sc.pp.highly_variable_genes(scranscanoramaadata, flavor='cell_ranger', n_top_genes=4000)
sc.pp.pca(scranscanoramaadata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(scranscanoramaadata, random_state = 2105, use_rep = "SC")
sc.tl.umap(scranscanoramaadata, random_state = 2105, n_components=3)

# # 5) Technical correction: Batch Correction using TrVAE
# rawadatas = adatas.copy()
# trvaecpmadata = 

# # 5.1) Normalizing & Extracting Top 1000 Highly Variable Genes
# # We can preserve more genes (up to 7000 like scGen) but in order to train the network quickly, we will extract top 1000 genes
# adata2 = adata.copy(); sc.pp.log1p(adata2); sc.pp.highly_variable_genes(adata2, n_top_genes=1000); adata2 = adata2[:, adata2.var['highly_variable']]
# # 5.2) Train/Test Split
# train_adata, valid_adata = trvae.utils.train_test_split(adata2, train_frac=0.80); train_adata.shape, valid_adata.shape # ((4348, 7109), (1087, 7109))
# # 5.3) Calculate number of batches
# n_conditions = len(train_adata.obs[condition_key].unique().tolist()); conditions = adata2.obs[condition_key].unique().tolist(); condition_encoder = trvae.utils.create_dictionary(conditions, []); condition_encoder; # {'bulk1001': 0, 'bulk997': 1, 'bulk1018': 2, 'stomach1001': 3}
# # 5.4)Create the network
# network = trvae.archs.trVAE(x_dimension=train_adata.shape[1], architecture=[128, 32], z_dimension=10, n_conditions=n_conditions, alpha=0.00005, beta=50, eta=100,  output_activation='relu')
# # 5.5) Training trVAE
# # NOTE: This will take 8 to 10 minutes on WS4 with single CPU
# outModelFile    = "{0}/trVAE_{1}_network.model".format(dataDir,projName)
# network.train(train_adata, valid_adata, condition_encoder, condition_key, n_epochs=500, batch_size=1024, verbose=5, early_stop_limit=300, lr_reducer=0, shuffle=True, save=True)
# # 5.6) Getting batch-corrected adata
# labels, _ = trvae.tl.label_encoder(adata2, condition_key=condition_key, label_encoder=condition_encoder); network.get_corrected(adata2, labels, return_z=True)
# # 5.7) MMD Layer UMAP visualization
# adata.obsm['mmd_latent']    = adata2.obsm['mmd_latent']; adata.obsm['z_latent']      = adata2.obsm['z_latent']; adata.obsm['reconstructed'] = adata2.obsm['reconstructed']
# sc.pp.neighbors(adata, random_state = 2105, n_neighbors=7, use_rep = "mmd_latent"); sc.tl.umap(adata, random_state = 2105, n_components=3); fig = plt.figure(figsize=(16,13)); fig.suptitle('TissueID')

# # 2D projection
# ax = fig.add_subplot(2, 2, 1);                  sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw UMAP")
# ax = fig.add_subplot(2, 2, 2);                  sc.pl.umap(adata, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE UMAP")
# # 3D projection
# ax = fig.add_subplot(2, 2, 3, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw UMAP")
# ax = fig.add_subplot(2, 2, 4, projection='3d'); sc.pl.umap(adata                , ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE UMAP")
# plt.tight_layout()
# plt.savefig("{0}/03_norm_TrVAE_batchCorrection_{1}_tissueID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

####
fig = plt.figure(figsize=(32,8))
# 2D projection
ax = fig.add_subplot(2, 5, 1);                  sc.pl.umap(adata             ,                  ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 5, 2);                  sc.pl.umap(cpmcombatadata   , legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="CPM Combat UMAP")
ax = fig.add_subplot(2, 5, 3);                  sc.pl.umap(scrancombatadata         , legend_loc=None, ax=ax, color="batch"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Scran Combat UMAP")
ax = fig.add_subplot(2, 5, 4);                  sc.pl.umap(cpmscanoramaadata, legend_loc=None, ax=ax, color="batch"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="CPM Scanorama UMAP")
ax = fig.add_subplot(2, 5, 5);                  sc.pl.umap(scranscanoramaadata, legend_loc=None, ax=ax, color="batch"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Scran Scanorama UMAP")
# 3D projection
ax = fig.add_subplot(2, 5, 6, projection='3d'); sc.pl.umap(adata             , legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 5, 7, projection='3d'); sc.pl.umap(cpmcombatadata   , legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="CPM Combat UMAP")
ax = fig.add_subplot(2, 5, 8, projection='3d'); sc.pl.umap(scrancombatadata         , legend_loc=None, ax=ax, color="batch"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Scran Combat UMAP")
ax = fig.add_subplot(2, 5, 9, projection='3d'); sc.pl.umap(cpmscanoramaadata, legend_loc=None, ax=ax, color="batch"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="CPM Scanorama UMAP")
ax = fig.add_subplot(2, 5, 10, projection='3d'); sc.pl.umap(scranscanoramaadata, legend_loc=None, ax=ax, color="batch"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Scran Scanorama UMAP")
plt.tight_layout()
plt.savefig("{0}/03_norm_all_batchCorrection_{1}_tissueID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')


# 4.4) Save the normalized cell cycle corrected adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/03_norm_all_batchCorrection_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/03_norm_all_batchCorrection_adata.h5ad" .format(dataDir, projName); normadata  = sc.read_h5ad(adatafile)
# normadata = adata.copy()

#########################################################################


# 5) Clustering
# 5.1) Perform clustering - using highly variable genes

############################################################
############################################################
# Using SCANORAMA for batch correction
############################################################
rawadatafile = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName)
rawadata     = sc.read_h5ad(rawadatafile)

# for adataLabel, normbcadata in zip(['cpm_combat', 'scran_combat', 'cpm_scanorama', 'scran_scanorama'], [cpmcombatadata, scrancombatadata, cpmscanoramaadata, scranscanoramaadata]):
for adataLabel, normbcadata in zip(['cpm_scanorama', 'scran_scanorama'], [cpmscanoramaadata, scranscanoramaadata]):
    adata = normbcadata.copy()
    bname = "{0}_{1}".format(projName, adataLabel)
    print(adata.shape)
    print(bname)
    print()

    # 7) Clustering
    # 7.1) Perform clustering - using highly variable genes
    sc.tl.louvain(adata, key_added='louvain', random_state=2105)
    sc.tl.louvain(adata, resolution=1, key_added='louvain_r1', random_state=2105)
    sc.tl.louvain(adata, resolution=1.5, key_added='louvain_r1.5', random_state=2105)
    sc.tl.louvain(adata, resolution=2.0, key_added='louvain_r2', random_state=2105)

    for i in np.linspace(0.1,0.9,9):
        try:
            sc.tl.louvain(adata, resolution=i, key_added='louvain_r{0}'.format(i), random_state=2105)
            print(adata.obs['louvain_r{0:0.1f}'.format(i)].value_counts())
        except:
            print("- Error in r: {0}".format(i))
    sc.tl.louvain(adata, resolution=0.3, key_added='louvain_r0.3', random_state=2105)
    sc.tl.louvain(adata, resolution=0.7, key_added='louvain_r0.7', random_state=2105)

    # 4.3) Visualizations

    # Calculations for the visualizations
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(adata, random_state = 2105, use_rep = "SC")
    sc.tl.umap(adata, random_state = 2105, n_components=3)

    # Plot visualizations
    # Visualize the clustering and how this is reflected by different technical covariates
    sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
    plt.savefig("{0}/03_{1}_clustering_all_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
    sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
    plt.savefig("{0}/03_{1}_clustering_all_louvain_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

    # Subcluster keys
    cluster_key            = "louvain_r0.5"
    cluster_bname          = "louvain_r05"
    subplot_title_fontsize = 12
    subplot_title_width    = 50

    # Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
    cluster_key_groups = adata.obs[cluster_key].cat.categories.tolist()
    cluster_cell_count = adata.obs[cluster_key].value_counts().to_dict()

    # Louvain UMAPs
    fig = plt.figure(figsize=(24,8))
    # 2D projection
    ax = fig.add_subplot(2, 4, 1);                  sc.pl.umap(adata, ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="louvain_r0.5 UMAP")
    ax = fig.add_subplot(2, 4, 2);                  sc.pl.umap(adata, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="tissueID UMAP")
    ax = fig.add_subplot(2, 4, 3);                  sc.pl.umap(adata, legend_loc=None, ax=ax, color="log_counts"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="log_counts UMAP")
    ax = fig.add_subplot(2, 4, 4);                  sc.pl.umap(adata, legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="mt_frac UMAP")
    # 3D projection
    ax = fig.add_subplot(2, 4, 5, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="louvain_r0.5 UMAP")
    ax = fig.add_subplot(2, 4, 6, projection='3d'); sc.pl.umap(adata, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="tissueID UMAP")
    ax = fig.add_subplot(2, 4, 7, projection='3d'); sc.pl.umap(adata         , legend_loc=None, ax=ax, color="log_counts"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="log_counts UMAP")
    ax = fig.add_subplot(2, 4, 8, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="mt_frac UMAP")
    plt.tight_layout()
    plt.savefig("{0}/03_{1}_louvain_tissueID_counts_mtfrac_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

    # UMAPS
    fig = plt.figure(figsize=(16,6))
    fig.suptitle(cluster_key)
    # 2D projection
    ax = fig.add_subplot(1, 2, 1);                  
    sc.pl.umap(adata, legend_loc=None, ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
    # 3D projection
    ax = fig.add_subplot(1, 2, 2, projection='3d'); 
    sc.pl.umap(adata, ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
    plt.savefig("{0}/03_{1}_clustering_{2}_UMAP_2D3D.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

    # 9.3) Plot separate bar plots, coloured in by cluster annotation, for each tissue
    # Convert palette into colormap
    clcmap = ListedColormap(sc.pl.palettes.vega_20)
    # Get the DF of tissue and clusters
    clusterBatchDF = adata.obs[['batch','{0}'.format(cluster_key)]].copy()
    # Replace batch number with batch names
    clusterBatchDF.replace({'batch': tissueIdDict}, inplace=True)
    # Remove index for groupby
    clusterBatchDF.reset_index(drop=True, inplace=True)
    # Get the number of cells for each cluster in every tissue
    ncellsClusterBatchDF = clusterBatchDF.groupby(['batch','{0}'.format(cluster_key)]).size()
    # Get the percent of cells for each cluster in every tissue 
    pcellsClusterBatchDF = pd.crosstab(index=clusterBatchDF['batch'], columns=clusterBatchDF['{0}'.format(cluster_key)], values=clusterBatchDF['{0}'.format(cluster_key)], aggfunc='count', normalize='index')
    # Plot the barplots
    fig = plt.figure(figsize=(16,6)); fig.suptitle("Cells for each {0} in each tissue".format(cluster_key))
    # plot numbers of cells
    ax = fig.add_subplot(1, 2, 1); ncellsClusterBatchDF.unstack().plot(kind='barh', stacked=True, colormap=clcmap, ax=ax, legend=None, title="Number of cells")
    # plot percent of cells
    ax = fig.add_subplot(1, 2, 2); pcellsClusterBatchDF.plot(kind='barh',stacked=True, colormap=clcmap, ax=ax, title="% of cells")
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='{0}'.format(cluster_key), title_fontsize=12)
    plt.tight_layout() # For non-overlaping subplots
    plt.savefig("{0}/03_{1}_clustering_{2}_tissueID_cluster_barplot.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')
