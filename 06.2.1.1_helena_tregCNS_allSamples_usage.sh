
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
# Sample quality plots
ax = fig.add_subplot(3, 2, 1); t1 = sc.pl.violin(adata, ['n_genes_by_counts', 'n_counts'], jitter=0.4, size=2, log=True, cut=0, ax = ax, show=False)
ax = fig.add_subplot(3, 2, 2); t2 = sc.pl.violin(adata, ['mt_frac','rb_frac'], jitter=0.4, size=2, log=False, cut=0, ax = ax, show=False)
# 1.2.4) Thresholdingecision based on counts
ax = fig.add_subplot(3, 2, 3); p3 = sns.distplot(adata.obs['n_counts'], kde=False, ax = ax, bins=50); #plt.show()
ax = fig.add_subplot(3, 2, 4); p4 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<2000], kde=False, ax = ax, bins=50); #plt.show()
# 1.2.5) Thresholding decision based on genes
ax = fig.add_subplot(3, 2, 5); p6 = sns.distplot(adata.obs['n_genes'], kde=False, ax = ax, bins=50); # plt.show()
ax = fig.add_subplot(3, 2, 6); p7 = sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, ax = ax, bins=50); # plt.show()
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

# 4.2) Batch Correction using Scanorama
scranadata2 = sc.AnnData(X=scranadata.X, var=scranadata.var, obs = scranadata.obs)
#variable genes for the full dataset
sc.pp.highly_variable_genes(scranadata2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'batch')
var_genes_batch = scranadata2.var.highly_variable_nbatches > 0
var_select = scranadata2.var.highly_variable_nbatches > 1
var_genes = var_select.index[var_select]
# Split per batch into new objects.
batches = ['0','1','2','3','4','5','6','7','8','9','10','11','12']
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
adatafile  = "{0}/02_normCC_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/02_normCC_{1}_adata.h5ad" .format(dataDir, projName); normadata  = sc.read_h5ad(adatafile)
# normadata = adata.copy()

#########################################################################
rawadatafile = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName)
rawadata     = sc.read_h5ad(rawadatafile)
# adata = normadata.copy()

# 5) Clustering
# 5.1) Perform clustering - using highly variable genes