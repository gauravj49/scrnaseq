# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

# Preprocessing
# Transpose the raw genesymbols file for genes
datamash -H transpose < /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/hippocampal/12weeks_hippo_publicDS3_counts_raw_genesymbols.txt > /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/hippocampal/T_12weeks_hippo_publicDS3_counts_raw_genesymbols.txt

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
projName        = "hippocampal" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/{0}/01_analysis1".format(projName); create_dir("{0}".format(output_dir))
cc_genes_file   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/Macosko_cell_cycle_genes.txt"
minGenesPerCell = 100
minCountPerCell = 300
minCellsPergene = 25
bname           = projName
plotsDir        = "{0}/plots".format(output_dir); create_dir(plotsDir)
dataDir         = "{0}/data".format(output_dir) ; create_dir(dataDir)

# 1 Reading the data
adata = sc.read_text('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/hippocampal/T_12weeks_hippo_publicDS3_counts_raw_genesymbols.txt')

# Make variable names unique
adata.var_names_make_unique()

# add celltype to the data
adata.obs['cellType'] = ['CA23' if 'CA23' in x else 'CA1' if 'CA1' in x else 'DG'  for x in adata.obs_names]

# Checking the total size of the data set
adata.shape # We have 521 cells and 48376 genes in the dataset

# Make a copy of the original data
origadata = adata.copy()

# AnnData object with n_obs × n_vars = 521 × 48376 
#     obs: 'batch', 'tissueID'
#     var: 'gene_ids', 'feature_types', 'genome'

# Checking the total size of the data set
adata.shape # We have 521 cells and 48376 genes in the dataset

# 2) Quality control 
# 2.1) Calculate QC covariates
adata.obs['n_counts']   = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes']    = (adata.X > 0).sum(1)
mt_gene_mask            = [gene.startswith('mt-') for gene in adata.var_names]
adata.obs['mt_frac']    = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']

# 2.2) Plot QC metrics
# Sample quality plots
t1 = sc.pl.violin(adata, 'n_counts', groupby='cellType', size=2, log=True, cut=0, show=False)
plt.savefig("{0}/01_raw_{1}_cellType_nCounts_plot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
t2 = sc.pl.violin(adata, 'mt_frac', groupby='cellType', show=False)
plt.savefig("{0}/01_raw_{1}_cellType_mtFraction_plot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.3) Data quality summary plots
p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/01_raw_{1}_genes_counts_scatterplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# p2 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', show=False)
# plt.savefig("{0}/01_{1}_genes_counts_scatterplot_zoomedin.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.3) Thresholding decision based on counts
p3 = sns.distplot(adata.obs['n_counts'], kde=False); #plt.show()
plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p4 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<2000], kde=False, bins=50); #plt.show()
plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot_lessthan_2000.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p5 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>5000], kde=False, bins=50); #plt.show()
plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot_greaterthan_5000.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.4) Thresholding decision based on genes
p6 = sns.distplot(adata.obs['n_genes'], kde=False, bins=50); # plt.show()
plt.savefig("{0}/01_raw_{1}_genes_histogramplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p7 = sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, bins=50); # plt.show()
plt.savefig("{0}/01_raw_{1}_genes_histogramplot_lessthan_1000.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.5) Filter cells according to identified QC thresholds:
origadata = adata.copy()
print('Total number of cells: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_counts = minCountPerCell)
print('Number of cells after min count filter: {:d}'.format(adata.n_obs))

# sc.pp.filter_cells(adata, max_counts = 40000)
# print('Number of cells after max count filter: {:d}'.format(adata.n_obs))

adata = adata[adata.obs['mt_frac'] < 0.05]
print('Number of cells after MT filter: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_genes = minGenesPerCell)
print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

# Total number of cells: 521
# Number of cells after min count filter: 521
# Number of cells after MT filter: 521
# Trying to set attribute `.obs` of view, copying.
# Number of cells after gene filter: 521

# 2.6) Filter genes according to identified QC thresholds:
# Min 5 cells - filters out 0 count genes
print('Total number of genes: {:d}'.format(adata.n_vars))
sc.pp.filter_genes(adata, min_cells=minCellsPergene)
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))

# Total number of genes: 55488
# Number of genes after cell filter: 14907

# 2.8) Calculations for the visualizations
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(adata, random_state = 2105)
sc.tl.umap(adata, random_state = 2105, n_components=3)

# 2.9) Plot visualizations
sc.pl.pca_scatter(adata, color='n_counts',show=False)
plt.savefig("{0}/01_raw_{1}_clustering_ncounts_PCA.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adata, color=['cellType'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/01_raw_{1}_clustering_cellType_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adata, color=['cellType'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/01_raw_{1}_clustering_cellType_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

########################
# 3) Normalization using scran
# Normalization using SCRAN
# This method requires a coarse clustering input to improve size factor esimation performance. 
# Thus, we use a simple preprocessing approach and cluster the data at a low resolution to get an 
# input for the size factor estimation. The basic preprocessing includes assuming all size factors 
# are equal (library size normalization to counts per million - CPM) and log-transforming the count data.

# 3.1) Perform a clustering for scran normalization in clusters
adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
sc.pp.neighbors(adata_pp)
sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)

# 3.2) Preprocess variables for scran normalization
input_groups = adata_pp.obs['groups']
data_mat = adata.X.T

# 3.3) Run scran in R
%%R -i data_mat -i input_groups -o size_factors
size_factors = computeSumFactors(data_mat, clusters=input_groups, min.mean=0.1)

# Delete adata_pp
del adata_pp

# 3.4) Visualize the estimated size factors
adata.obs['size_factors'] = size_factors
sc.pl.scatter(adata, 'size_factors', 'n_counts', show=False)
plt.savefig("{0}/02_norm_{1}_sizefactors_vs_ncounts.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.scatter(adata, 'size_factors', 'n_genes' , show=False)
plt.savefig("{0}/02_norm_{1}_sizefactors_vs_ngenes.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sns.distplot(size_factors, bins=50, kde=False)
plt.savefig("{0}/02_norm_{1}_sizefactors_histogram.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Keep the count data in a counts layer
adata.layers["counts"] = adata.X.copy()

# 3.5) Normalize adata 
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)

# Store the full data set in 'raw' as log-normalised data for statistical testing
adata.raw = adata

# Highly variable gene information is stored automatically in the adata.var['highly_variable'] field. The dataset now contains:
# a 'counts' layer with count data
# log-normalized data in adata.raw
# batch corrected data in adata.X
# highly variable gene annotations in adata.var['highly_variable']
# The HVG labels will be used to subselect genes for clustering and trajectory analysis.
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))

# 3.6) Biological correction
# Score cell cycle and visualize the effect:
cc_genes         = pd.read_table(cc_genes_file, delimiter='\t')
s_genes          = cc_genes['S'].dropna()
g2m_genes        = cc_genes['G2.M'].dropna()
s_genes_mm       = [gene.lower().capitalize() for gene in s_genes]
g2m_genes_mm     = [gene.lower().capitalize() for gene in g2m_genes]
s_genes_mm_ens   = adata.var_names[np.in1d(adata.var_names, s_genes_mm)]
g2m_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, g2m_genes_mm)]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)

# 3.7) Calculations for the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(adata, random_state = 2105)
sc.tl.umap(adata, random_state = 2105, n_components=3)

# 3.8) Plot visualizations
sc.pl.umap(adata, color=['S_score', 'G2M_score'], use_raw=False, palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_S_G2M_Phase_UMAP.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.umap(adata, color='phase', use_raw=False, palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_combined_Phase_UMAP.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.pca_scatter(adata, color='n_counts',show=False)
plt.savefig("{0}/02_norm_{1}_clustering_ncounts_PCA.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.umap(adata, color=['cellType'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_cellType_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.umap(adata, color=['cellType'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/02_norm_{1}_clustering_cellType_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')


# 5) Clustering
# 5.1) Perform clustering - using highly variable genes
sc.tl.louvain(adata, key_added='louvain'     , random_state=2105)
sc.tl.louvain(adata, key_added='louvain_r1'  , random_state=2105, resolution=1.0)
sc.tl.louvain(adata, key_added='louvain_r1.5', random_state=2105, resolution=1.5)
sc.tl.louvain(adata, key_added='louvain_r2'  , random_state=2105, resolution=2.0)

for i in np.linspace(0.1,0.9,9):
  i = round(i,2)    
  try:
      sc.tl.louvain(adata, resolution=i, key_added='louvain_r{0}'.format(i), random_state=2105)
      print(adata.obs['louvain_r{0:0.1f}'.format(i)].value_counts())
      print()
  except:
      print("- Error in r: {0}".format(i))

# 0    323
# 1    159
# 2     39
# Name: louvain_r0.1, dtype: int64

# 0    322
# 1    160
# 2     39
# Name: louvain_r0.2, dtype: int64

# 0    289
# 1    161
# 2     39
# 3     32
# Name: louvain_r0.3, dtype: int64

# 0    289
# 1    101
# 2     60
# 3     39
# 4     32
# Name: louvain_r0.4, dtype: int64

# 0    290
# 1    100
# 2     60
# 3     39
# 4     32
# Name: louvain_r0.5, dtype: int64

# 0    290
# 1     95
# 2     65
# 3     39
# 4     32
# Name: louvain_r0.6, dtype: int64

# 0    202
# 1     95
# 2     88
# 3     65
# 4     39
# 5     32
# Name: louvain_r0.7, dtype: int64

# 0    202
# 1     95
# 2     88
# 3     65
# 4     39
# 5     32
# Name: louvain_r0.8, dtype: int64

# 0    108
# 1     95
# 2     93
# 3     88
# 4     65
# 5     39
# 6     33
# Name: louvain_r0.9, dtype: int64


print(adata.obs['louvain'].value_counts())
print(adata.obs['louvain_r1'].value_counts())
print(adata.obs['louvain_r1.5'].value_counts())
print(adata.obs['louvain_r2'].value_counts())

# 0    95
# 1    93
# 2    89
# 3    83
# 4    65
# 5    39
# 6    32
# 7    25
# Name: louvain, dtype: int64
# 0    95
# 1    93
# 2    89
# 3    83
# 4    65
# 5    39
# 6    32
# 7    25
# Name: louvain_r1, dtype: int64
# 0    96
# 1    92
# 2    88
# 3    69
# 4    67
# 5    42
# 6    39
# 7    28
# Name: louvain_r1.5, dtype: int64
# 0     93
# 1     66
# 2     55
# 3     54
# 4     47
# 5     41
# 6     39
# 7     36
# 8     32
# 9     31
# 10    23
# 11     4
# Name: louvain_r2, dtype: int64

# Plot visualizations
# Visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=['louvain'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_all_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Looking at the data, choose r = 0.5
sc.pl.umap(adata, color=['louvain_r0.5'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_louvain_r05_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.umap(adata, color=['louvain_r0.5'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, legend_loc='on data', show=False)
plt.savefig("{0}/02_norm_{1}_clustering_louvain_r05_legend_onData_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
























# 6) Clustering
# 7.1) Perform clustering - using highly variable genes
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
sc.tl.louvain(adata, key_added='louvain', random_state=2105)

for i in np.linspace(0.1,0.9,9):
  print(adata.obs['louvain_r{0:0.1f}'.format(i)].value_counts())


# Number of cells in each cluster
adata.obs['louvain_r0.5'].value_counts()
# 0    2232
# 1    1197
# 2    1090
# 3     864
# 4     831
# 5     504
# 6     226
# 7     207
# 8     110
# Name: louvain_r0.5, dtype: int64

# 4.3) Visualizations

# Calculations for the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata, random_state = 2105)
# sc.tl.diffmap(adata)
# sc.tl.draw_graph(adata)

# Plot visualizations
# Visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=['louvain', 'louvain_r0.5'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/{1}_clustering_louvain_r05_UMAP.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')

# UMAPs
sc.pl.pca_scatter(adata, color='n_counts', palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_ncounts_PCA.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['tissueID'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_tissueID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color='n_counts', palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_ncounts_UMAP.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['log_counts', 'mt_frac'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_logCounts_mtFrac_UMAP.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')

# sc.pl.pca_scatter(adata, color='n_counts')
# sc.pl.umap(adata, color='n_counts')
# sc.pl.umap(adata, color='tissueID')
# sc.pl.diffmap(adata, color='n_counts', components=['1,2','1,3'])
# sc.pl.draw_graph(adata, color='n_counts')

# At a resolution of r=0.5 the broad clusters in the visualization are captured well in the data. 
# The covariate plots show that clusters 0 and 6 in this data set are characterized by low and high 
# counts respectively. In the case of cluster 6 this may be biologically relevant, while cluster 0 
# is also characterized by higher mitochondrial read fractions. This indicates cell stress.

# 7.2) Marker genes & cluster annotation
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby='louvain_r0.5', key_added='rank_genes_r0.5')
sc.tl.rank_genes_groups(adata, groupby='louvain_r1', key_added='rank_genes_r1')

# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_r0.5', groups=['0','1','2'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_0_1_2.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.rank_genes_groups(adata, key='rank_genes_r0.5', groups=['3','4','5'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_3_4_5.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.rank_genes_groups(adata, key='rank_genes_r0.5', groups=['6', '7', '8'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_7_6_8.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')

# Pairwise between 0,2,5,10
for l in list(itertools.combinations(['0','2','5','10'], 2)):
  grp = l[0]
  ref = l[1]
  print(list(l))
  sc.tl.rank_genes_groups(adata, groupby='louvain_r1', key_added='rank_genes_{0}_{1}'.format(grp, ref), groups = list(grp), reference = '{}'.format(ref))
  sc.pl.rank_genes_groups(adata, key='rank_genes_{0}_{1}'.format(grp, ref), fontsize=12, show=False)
  plt.savefig("{0}/{1}_marker_genes_ranking_pairwise_{2}_{3}.png".format(plotsDir, bname, grp, ref) , bbox_inches='tight'); plt.close('all')


# Known marker genes:
# marker_genes = dict()
# marker_genes['Stem'] = ['Lgr5', 'Ascl2', 'Slc12a2', 'Axin2', 'Olfm4', 'Gkn3']
# marker_genes['Enterocyte (Proximal)'] = ['Gsta1','Rbp2','Adh6a','Apoa4','Reg3a','Creb3l3','Cyp3a13','Cyp2d26','Ms4a10','Ace','Aldh1a1','Rdh7','H2-Q2', 'Hsd17b6','Gstm3','Gda','Apoc3','Gpd1','Fabp1','Slc5a1','Mme','Cox7a1','Gsta4','Lct','Khk','Mttp','Xdh','Sult1b1', 'Treh','Lpgat1','Dhrs1','Cyp2c66','Ephx2','Cyp2c65','Cyp3a25','Slc2a2','Ugdh','Gstm6','Retsat','Ppap2a','Acsl5', 'Cyb5r3','Cyb5b','Ckmt1','Aldob','Ckb','Scp2','Prap1']
# marker_genes['Enterocyte (Distal)'] = ['Tmigd1','Fabp6','Slc51b','Slc51a','Mep1a','Fam151a','Naaladl1','Slc34a2','Plb1','Nudt4','Dpep1','Pmp22','Xpnpep2','Muc3','Neu1','Clec2h','Phgr1','2200002D01Rik','Prss30','Cubn','Plec','Fgf15','Crip1','Krt20','Dhcr24','Myo15b','Amn','Enpep','Anpep','Slc7a9','Ocm','Anxa2','Aoc1','Ceacam20','Arf6','Abcb1a','Xpnpep1','Vnn1','Cndp2','Nostrin','Slc13a1','Aspa','Maf','Myh14']
# marker_genes['Goblet'] = ['Agr2', 'Fcgbp', 'Tff3', 'Clca1', 'Zg16', 'Tpsg1', 'Muc2', 'Galnt12', 'Atoh1', 'Rep15', 'S100a6', 'Pdia5', 'Klk1', 'Pla2g10', 'Spdef', 'Lrrc26', 'Ccl9', 'Bace2', 'Bcas1', 'Slc12a8', 'Smim14', 'Tspan13', 'Txndc5', 'Creb3l4', 'C1galt1c1', 'Creb3l1', 'Qsox1', 'Guca2a', 'Scin', 'Ern2', 'AW112010', 'Fkbp11', 'Capn9', 'Stard3nl', 'Slc50a1', 'Sdf2l1', 'Hgfa', 'Galnt7', 'Hpd', 'Ttc39a', 'Tmed3', 'Pdia6', 'Uap1', 'Gcnt3', 'Tnfaip8', 'Dnajc10', 'Ergic1', 'Tsta3', 'Kdelr3', 'Foxa3', 'Tpd52', 'Tmed9', 'Spink4', 'Nans', 'Cmtm7', 'Creld2', 'Tm9sf3', 'Wars', 'Smim6', 'Manf', 'Oit1', 'Tram1', 'Kdelr2', 'Xbp1', 'Serp1', 'Vimp', 'Guk1', 'Sh3bgrl3', 'Cmpk1', 'Tmsb10', 'Dap', 'Ostc', 'Ssr4', 'Sec61b', 'Pdia3', 'Gale', 'Klf4', 'Krtcap2', 'Arf4', 'Sep15', 'Ssr2', 'Ramp1', 'Calr', 'Ddost']
# marker_genes['Paneth'] = ['Gm15284', 'AY761184', 'Defa17', 'Gm14851', 'Defa22', 'Defa-rs1', 'Defa3', 'Defa24', 'Defa26', 'Defa21', 'Lyz1', 'Gm15292', 'Mptx2', 'Ang4']
# marker_genes['Enteroendocrine'] = ['Chgb', 'Gfra3', 'Cck', 'Vwa5b2', 'Neurod1', 'Fev', 'Aplp1', 'Scgn', 'Neurog3', 'Resp18', 'Trp53i11', 'Bex2', 'Rph3al', 'Scg5', 'Pcsk1', 'Isl1', 'Maged1', 'Fabp5', 'Celf3', 'Pcsk1n', 'Fam183b', 'Prnp', 'Tac1', 'Gpx3', 'Cplx2', 'Nkx2-2', 'Olfm1', 'Vim', 'Rimbp2', 'Anxa6', 'Scg3', 'Ngfrap1', 'Insm1', 'Gng4', 'Pax6', 'Cnot6l', 'Cacna2d1', 'Tox3', 'Slc39a2', 'Riiad1']
# marker_genes['Tuft'] = ['Alox5ap', 'Lrmp', 'Hck', 'Avil', 'Rgs13', 'Ltc4s', 'Trpm5', 'Dclk1', 'Spib', 'Fyb', 'Ptpn6', 'Matk', 'Snrnp25', 'Sh2d7', 'Ly6g6f', 'Kctd12', '1810046K07Rik', 'Hpgds', 'Tuba1a', 'Pik3r5', 'Vav1', 'Tspan6', 'Skap2', 'Pygl', 'Ccdc109b', 'Ccdc28b', 'Plcg2', 'Ly6g6d', 'Alox5', 'Pou2f3', 'Gng13', 'Bmx', 'Ptpn18', 'Nebl', 'Limd2', 'Pea15a', 'Tmem176a', 'Smpx', 'Itpr2', 'Il13ra1', 'Siglecf', 'Ffar3', 'Rac2', 'Hmx2', 'Bpgm', 'Inpp5j', 'Ptgs1', 'Aldh2', 'Pik3cg', 'Cd24a', 'Ethe1', 'Inpp5d', 'Krt23', 'Gprc5c', 'Reep5', 'Csk', 'Bcl2l14', 'Tmem141', 'Coprs', 'Tmem176b', '1110007C09Rik', 'Ildr1', 'Galk1', 'Zfp428', 'Rgs2', 'Inpp5b', 'Gnai2', 'Pla2g4a', 'Acot7', 'Rbm38', 'Gga2', 'Myo1b', 'Adh1', 'Bub3', 'Sec14l1', 'Asah1', 'Ppp3ca', 'Agt', 'Gimap1', 'Krt18', 'Pim3', '2210016L21Rik', 'Tmem9', 'Lima1', 'Fam221a', 'Nt5c3', 'Atp2a3', 'Mlip', 'Vdac3', 'Ccdc23', 'Tmem45b', 'Cd47', 'Lect2', 'Pla2g16', 'Mocs2', 'Arpc5', 'Ndufaf3']

# Read the marker genes into a pandas dataframe
marker_file  = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_V1.txt'
markersDF    = pd.read_csv(marker_file, sep="\t")
marker_genes = markersDF.groupby('CellLines')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()
marker_genes_cellTypes = markersDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()

# For mouse cell atlas marker genes
ma_marker_file       = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_mouse_cellatlas_V1.txt'
ma_markersDF         = pd.read_csv(ma_marker_file, sep="\t", header=None, index_col=None)
ma_markersDF         = ma_markersDF[0].str.split(",", n = 1, expand = True)
ma_markersDF.columns = ['CellTypes', 'MarkerGenes']
ma_marker_genes      = ma_markersDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()


cell_annotation = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_r0.5')
ma_cell_annotation = sc.tl.marker_gene_overlap(adata, ma_marker_genes, key='rank_genes_r0.5')
cell_annotation
#                  0    1    2    3    4    5    6    7    8    9   10   11
# Endocrine      0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
# Epithelial     0.0  0.0  0.0  0.0  0.0  1.0  0.0  1.0  0.0  1.0  0.0  0.0
# Immune         0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  2.0  0.0  0.0  1.0
# Intestinal     0.0  0.0  0.0  2.0  4.0  0.0  1.0  0.0  0.0  0.0  0.0  1.0
# Liver          0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
# Mucous         0.0  0.0  0.0  1.0  1.0  0.0  2.0  0.0  0.0  0.0  0.0  0.0
# Myocytes       0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
# Pancreas       0.0  0.0  4.0  5.0  2.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
# Proliferative  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
# Stem           0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
# Stomach        0.0  0.0  0.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0
# Stroma         0.0  0.0  0.0  0.0  0.0  0.0  5.0  0.0  0.0  0.0  0.0  0.0

cell_annotation_norm = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_r0.5', normalize='reference')
sns.heatmap(cell_annotation_norm, cbar=False, annot=True)
plt.savefig("{0}/{1}_marker_genes_cell_annotation_norm_heatmap.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')

ma_cell_annotation_norm = sc.tl.marker_gene_overlap(adata, ma_marker_genes, key='rank_genes_r0.5', normalize='reference')
sns.heatmap(ma_cell_annotation_norm, cbar=False, annot=True)
plt.savefig("{0}/{1}_stomach_marker_list_mouse_cellatlas_V1_annotation_norm_heatmap.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')

# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

sc.pl.umap(adata, color=['louvain', 'Vim', 'Mdm2', 'Trp53', 'Irf8', 'Myc', 'Gamt', 'E2f1', 'Pcna', 'Tgfbr2'], use_raw=False, color_map=mymap)
plt.savefig("{0}/{1}_marker_genes_adult_stomach_mouse_cell_atlas_UMAPs.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')


sc.pl.umap(adata, color=['Ctrb1','Cd2', 'Myc', 'Rb1','Cdkn2a'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)
sc.pl.umap(adata, color=['Birc5', 'Casp3', 'Stat3', 'Alb'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

sc.pl.umap(adata, color=['Birc5', 'Stat3', 'Reg1', 'Gm26917', 'Ctrb1', 'Clps', 'Hbb-bs', 'Hba-a1','Hba-a2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)
sc.pl.umap(adata, color=['Birc5', 'Blvrb', 'Car2', 'Hbb-bt', 'Clps', 'Hbb-bs', 'Hba-a1','Hba-a2', 'Hspa1b', 'Apoe', 'C1qb', 'Cxcl2', 'Slpi'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)
sc.pl.umap(adata, color=['Sh2d1a', 'Cd3d','Cd3e','Cd8a','Retnlg','S100a8','S100a9','Cxcl2', 'Slpi', 'Srgn', 'Cd84', 'Stip1','Cd44', 'Jak1'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

sc.pl.umap(adata, color=['Btg1', 'Morf4l2','Marcks','Ptprc','BC005537','Ctnnb1','Ptma','AW112010', 'Hnrnpf', 'Hspa1b', 'Hnrnph1', 'Dazap2','Laptm5', 'Id2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

sc.pl.umap(adata, color=['louvain_r0.5','Apoe', 'Dcn', 'Cd74','Pf4', 'Lyz2', 'Ly6c1', 'Plvap', 'Fabp5', 'Birc5', 'Ube2c', 'Dmbt1', 'Cyr61', 'Igfbp3', 'Clps'], use_raw=False, color_map='hot', show=False)
plt.savefig("{0}/{1}_marker_genes_adult_stomach_mouse_cell_atlas_UMAPs.png".format(plotsDir, bname) , bbox_inches='tight'); plt.close('all')


# # Check expression of enterocyte markers
# #Collate all enterocyte markers and get the gene IDs in the data set
# ids_pancreas = np.in1d(adata.var_names, marker_genes['Pancreas'])

# #Calculate the mean expression of enterocyte markers
# adata.obs['Pancreas_marker_expr'] = adata.X[:,ids_pancreas].mean(1)

# #Plot enterocyte expression
# sc.pl.violin(adata, 'Pancreas_marker_expr', groupby='louvain_r0.5')
# sc.pl.umap(adata, color='Pancreas_marker_expr', color_map=mymap)

# Generate the UMAPs for each marker categories
for k in marker_genes.keys():
  ids = np.in1d(adata.var_names, marker_genes[k])
  adata.obs['{0}_marker_expr'.format(k)] = adata.X[:,ids].mean(1)
  sc.pl.violin(adata, '{0}_marker_expr'.format(k), groupby='louvain_r0.5', show=False)
  plt.savefig("{0}/{1}_marker_genes_stomach_{2}_violinPlots.png".format(plotsDir, bname, k) , bbox_inches='tight'); plt.close('all')
  sc.pl.umap(adata, color=['louvain','{0}_marker_expr'.format(k)], color_map=mymap, show=False)
  plt.savefig("{0}/{1}_marker_genes_stomach_{2}_UMAPs.png".format(plotsDir, bname, k) , bbox_inches='tight'); plt.close('all')
  

# Generate the UMAPs for each marker categories
plt.figure(figsize=(100,8))
for k in ma_marker_genes.keys():
  ids = np.in1d(adata.var_names, ma_marker_genes[k])
  adata.obs['{0}_ma_marker_expr'.format(k)] = adata.X[:,ids].mean(1)
  # sc.pl.violin(adata, '{0}_ma_marker_expr'.format(k), groupby='louvain_r0.5', show=False)
  # plt.savefig("{0}/{1}_mouse_cellatlas_marker_genes_stomach_{2}_violinPlots.png".format(plotsDir, bname, k) , bbox_inches='tight'); plt.close('all')
  # sc.pl.umap(adata, color=['louvain','{0}_ma_marker_expr'.format(k)], color_map=mymap, show=False)
  sc.pl.umap(adata, color=['{0}_ma_marker_expr'.format(k)], color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  plt.savefig("{0}/{1}_mouse_cellatlas_marker_genes_stomach_{2}_UMAPs.png".format(plotsDir, bname, k) , bbox_inches='tight', dpi=300); plt.close('all')
  
# Categories to rename
adata.obs['louvain_r0.5'].cat.categories

# On the subclusters
sbadata =adata.copy()
# sc.pl.umap(sbadata, color='louvain', palette=sc.pl.palettes.vega_20, size=50)
plt.figure(figsize=(100,8))
sc.pl.umap(sbadata, color='louvain_r1', palette=sc.pl.palettes.vega_20, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/{1}_Louvainr05_UMAPs.png".format(plotsDir, bname, k) , bbox_inches='tight', dpi=300); plt.close('all')

# Groups interested = 0, 2, 5, 10
sc.tl.rank_genes_groups(sbadata, groupby='louvain_r1', key_added='rank_genes_0_2', groups = '2', reference = '0')
sc.pl.rank_genes_groups(sbadata, key='rank_genes_0_2', fontsize=12)


grp = ['4', '6', '13']
ref = '11'
sc.tl.rank_genes_groups(adata, groupby='louvain_r1', key_added='rank_genes_{0}_{1}'.format(grp, ref), groups = list(grp), reference =  '{}'.format(ref))
sc.pl.rank_genes_groups(adata, key='rank_genes_{0}_{1}'.format(grp, ref), fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_pairwise_{2}_{3}.png".format(plotsDir, bname, grp, ref) , bbox_inches='tight'); plt.close('all')

# Plot umap for the top genes for each cluster
clusterDir = "{0}/clusterDir".format(plotsDir); create_dir(clusterDir)
for cid in [0, 2, 5, 10]:
  cluster = adata[(adata.obs['louvain_r1'] == '{0}'.format(cid))].to_df() 
  topGenes = np.mean(cluster, axis=0).sort_values(ascending=False)[0:25].index.tolist()
  print("\n- Cluster{0}: {1}".format(cid, topGenes))
  plt.figure(figsize=(100,8))
  sc.pl.umap(adata, color=topGenes, use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  plt.savefig("{0}/{1}_top_25_genes_cluster{2}_UMAPs.png".format(clusterDir, bname, cid) , bbox_inches='tight'); plt.close('all')
  
# Plot umap for the top marker genes for each cluster
markerDir = "{0}/markerDir/pdf".format(plotsDir); create_dir(markerDir)
for cid in range(len(adata.uns['rank_genes_r1']['names'])):
  topGenes = adata.uns['rank_genes_r1']['names'][cid].tolist()
  print("\n- Group{0}: {1}".format(cid, topGenes))
  # plt.figure(figsize=(100,8))
  # sc.pl.umap(adata, color=topGenes, use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  # plt.savefig("{0}/{1}_marker_genes_group{2}_UMAPs.png".format(markerDir, bname, cid) , bbox_inches='tight'); plt.close('all')
  sc.pl.umap(adata, color=topGenes, use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor=None, show=False)
  plt.savefig("{0}/{1}_marker_genes_group{2}_UMAPs.pdf".format(markerDir, bname, cid) , bbox_inches='tight'); plt.close('all')
  
# Assign the celltypes to clusters
adata.rename_categories('louvain_r1', ['Birc5⁺/basal', 'Acinar', 'Gastric mucosa/Basal', 'Tcells', 'Dendritic', '5', 'Macrophages', 'Endothelial/Epithelial_Igfbp3⁺', 'Erythrocytes', 'Fibroblasts', '10', 'Restin', ' ', 'Dendritic/Macrophages'])
adata.obs['louvain_r1'].value_counts()
plt.figure(figsize=(100,8))
sc.pl.umap(adata, color='louvain_r1', palette=sc.pl.palettes.vega_20, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/{1}_Clusters_CellTypes_UMAPs.png".format(plotsDir, bname, k) , bbox_inches='tight', dpi=300); plt.close('all')


# DO Not Delete ... this can be use to visualize a subset of clusters  
# subadata02 = sbadata[(sbadata.obs['louvain'] == '0') | (sbadata.obs['louvain'] == '2')]
# subadata02.shape
# sc.tl.louvain(subadata02, restrict_to=('louvain', ['0', '2']), resolution=0.2, key_added='louvain_0_2')
# sc.pl.umap(subadata02, color='louvain_0_2', palette=sc.pl.palettes.vega_20, size=50)



# 7) Feature selection
# 8) Dimensionality reduction


sc.tl.paga(adata, groups='louvain_r1')
sc.pl.paga_compare(adata)
sc.pl.paga(adata)

fig1, ax1 = plt.subplots()
sc.pl.umap(adata, size=40, ax=ax1, show=False)
sc.pl.paga(adata, pos=adata.uns['paga']['pos'], show=False, node_size_scale=10, node_size_power=1, ax=ax1, text_kwds={'alpha':0})
#plt.savefig('./figures/umap_paga_overlay_gut.pdf', dpi=300, format='pdf')
plt.show()

#########################################
# Save session
import dill
filename = "{0}/01_{1}_analysis_1.pkl".format(output_dir, projName)
dill.dump_session(filename)

# and to load the session again:
import dill
filename = "{0}/01_{1}_analysis_1.pkl".format(output_dir, projName)
dill.load_session(filename)

#########################################

#########################################
# Saving session with pickle was giving the error
# Write the adata and cadata object to file
adatafile  = "{0}/{1}_adata.h5ad" .format(output_dir, projName); adata.write(adatafile)

# Read back the corrected adata object
adatafile  = "{0}/{1}_adata.h5ad" .format(output_dir, projName); adata  = sc.read_h5ad(adatafile)
#########################################

