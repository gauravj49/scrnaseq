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
projName        = "klegerDuctalDiff_gse84133" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/02_gse84133"; create_dir("{0}".format(output_dir))
minGenesPerCell = 100
minCellsPergene = 2
bname           = projName
qcDir           = "{0}/analysis".format(output_dir); create_dir(qcDir)
countsDir       = "{0}/counts".format(output_dir); create_dir(countsDir)
input_matrix_file = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/02_gse84133/scMatrix.txt'

# Add 
# 1) Reading the data
# head -3 output/klegerDuctalDiff/02_gse84133/scMatrix.txt | column -t | perl -lane 'print "@F[-6..-1]"'
# ┌───────┬──────┬────┬─────────────────────────────┬─────────┬───────┐
# │ ZZEF1 │ ZZZ3 │ pk │           CellID            │ Cluster │ Batch │
# ├───────┼──────┼────┼─────────────────────────────┼─────────┼───────┤
# │     0 │    0 │  1 │ human1_lib1.final_cell_0001 │ acinar  │     1 │
# │     0 │    1 │  0 │ human1_lib1.final_cell_0002 │ acinar  │     1 │
# └───────┴──────┴────┴─────────────────────────────┴─────────┴───────┘
orgdataDF   = pd.read_csv(input_matrix_file, sep="\t", index_col=['CellID'])

orgcountsDF = orgdataDF.copy()
orgcountsDF.drop(['Cluster', 'Batch'], axis=1, inplace=True)
adata = sc.AnnData(orgcountsDF)

# Make variable names unique
adata.var_names_make_unique()

# Add annotation
adata.obs['Batch']   = orgdataDF['Batch'].astype(str)
adata.obs['Cluster'] = orgdataDF['Cluster'].astype(str)
# AnnData object with n_obs × n_vars = 8569 × 20125 
#     obs: 'Batch', 'Cluster'

# Checking the total size of the data set
adata.shape # We have 8569 cells and 20125 genes in the dataset

# 2) Quality control 
# 2.1) Calculate QC covariates
adata.obs['n_counts']   = adata.X.sum(1)
adata.obs['n_genes']    = (adata.X > 0).sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])

# 2.2) Plot QC metrics
# Sample quality plots
t1 = sc.pl.violin(adata, 'n_counts', groupby='Batch', size=2, log=True, cut=0, show=False)
plt.savefig("{0}/{1}_Batch_nCounts_plot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

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

# 2.5) Filter cells according to identified QC thresholds:
origadata = adata.copy()
print('Total number of cells: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_counts = 300)
print('Number of cells after min count filter: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_genes = 300)
print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

# Total number of cells: 8569
# Number of cells after min count filter: 8569
# Number of cells after gene filter: 8569

# 2.6) Filter genes according to identified QC thresholds:
# Min 5 cells - filters out 0 count genes
print('Total number of genes: {:d}'.format(adata.n_vars))
sc.pp.filter_genes(adata, min_cells=10)
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
# Total number of genes: 20125
# Number of genes after cell filter: 15117

# # 2.7) Save the filtered raw data as tab separated matrix 
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
plt.savefig("{0}/01_raw_{1}_ncounts_PCA.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['Batch'], show=False)
plt.savefig("{0}/01_raw_{1}_Batch_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['Cluster'], show=False)
plt.savefig("{0}/01_raw_{1}_Cluster_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

########################
# 3) Expression recovery (denoising), Normalization and log transformation
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

# # Delete adata_pp
# del adata_pp

# Visualize the estimated size factors
adata.obs['size_factors'] = size_factors
sc.pl.scatter(adata, 'size_factors', 'n_counts', show=False)
plt.savefig("{0}/{1}_scrna_sizefactors_vs_ncounts.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.scatter(adata, 'size_factors', 'n_genes' , show=False)
plt.savefig("{0}/{1}_scrna_sizefactors_vs_ngenes.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sns.distplot(size_factors, bins=50, kde=False)
plt.savefig("{0}/{1}_scrna_sizefactors_histogram.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Keep the count data in a counts layer
adata.layers["counts"] = adata.X.copy()

# Normalize adata 
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)

# Store the full data set in 'raw' as log-normalised data for statistical testing
adata.raw = adata

# 4) Technical correction
# 4.1) Batch Correction using ComBat
sc.pp.combat(adata, key='Batch')

# 4.2) Highly Variable Genes
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig("{0}/{1}_highly_variable_genes.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Save the normalized, log transformed, batch and cell cycle corrected data
normCorrectedDF            = adata.to_df()
normCorrectedDF['Batch']   = orgdataDF['Batch'].astype(str)
normCorrectedDF['Cluster'] = orgdataDF['Cluster']
normCorrectedDF.to_csv("{0}/02_normalizedRaw_T_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="PatId")
normCorrectedDF.T.to_csv("{0}/02_normalizedRaw_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol")

# 7) Clustering
# # 7.1) Perform clustering - using highly variable genes
# sc.tl.louvain(adata, resolution=1.0, key_added='louvain_r1'  , random_state=2105)
# sc.tl.louvain(adata, resolution=0.9, key_added='louvain_r0.9', random_state=2105)
# sc.tl.louvain(adata, resolution=0.8, key_added='louvain_r0.8', random_state=2105)
# sc.tl.louvain(adata, resolution=0.7, key_added='louvain_r0.7', random_state=2105)
# sc.tl.louvain(adata, resolution=0.6, key_added='louvain_r0.6', random_state=2105)
# sc.tl.louvain(adata, resolution=0.5, key_added='louvain_r0.5', random_state=2105)
# sc.tl.louvain(adata, resolution=0.4, key_added='louvain_r0.4', random_state=2105)
# sc.tl.louvain(adata, resolution=0.3, key_added='louvain_r0.3', random_state=2105)
# sc.tl.louvain(adata, resolution=0.2, key_added='louvain_r0.2', random_state=2105)
# sc.tl.louvain(adata, resolution=0.1, key_added='louvain_r0.1', random_state=2105)
# sc.tl.louvain(adata, key_added='louvain', random_state=2105)

# for i in np.linspace(0.1,0.9,9):
#   print(adata.obs['louvain_r{0:0.1f}'.format(i)].value_counts())


# # Number of cells in each cluster
# adata.obs['louvain'].value_counts()
# # 0    2232
# # 1    1197
# # 2    1090
# # 3     864
# # 4     831
# # 5     504
# # 6     226
# # 7     207
# # 8     110
# # Name: louvain_r0.5, dtype: int64

# 4.3) Visualizations

# Calculations for the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata, random_state = 2105)

# UMAPs
sc.pl.pca_scatter(adata, color='n_counts', show=False)
plt.savefig("{0}/02_norm_{1}_ncounts_PCA.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['Batch'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_Batch_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adata, color=['Cluster'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_Cluster_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adata, color='n_counts', palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_ncounts_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# 7.2) Marker genes & cluster annotation
# Calculate marker genes
# sc.tl.rank_genes_groups(adata, groupby='louvain_r0.5', key_added='rank_genes_r0.5')
# sc.tl.rank_genes_groups(adata, groupby='louvain_r1', key_added='rank_genes_r1')
# sc.tl.rank_genes_groups(adata, groupby='louvain', key_added='rank_genes')
sc.tl.rank_genes_groups(adata, groupby='Cluster', key_added='rank_genes_Clusters', n_genes=adata.shape[1])

# adata.obs['Cluster'].value_counts()                                                                                                                           Out[67]: 
# beta                  2525
# alpha                 2326
# ductal                1077
# acinar                 958
# delta                  601
# activated_stellate     284
# gamma                  255
# endothelial            252
# quiescent_stellate     173
# macrophage              55
# mast                    25
# epsilon                 18
# schwann                 13
# t_cell                   7
# Name: Cluster, dtype: int64

# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_Clusters', groups=['beta','alpha','ductal','acinar','delta','activated_stellate','gamma','endothelial','quiescent_stellate','macrophage','mast','epsilon','schwann','t_cell'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_Clusters.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

# Save UMAP of marker genes
sc.pl.umap(adata, color=['Cluster', 'KRT19', 'KRT7', 'CFTR', 'MUC1','BHLHA15', 'PTF1A', 'CPA1', 'AMY2A'], use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/{1}_marker_genes_ductal_acinar_UMAPs.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
for g in ['Cluster', 'KRT19', 'KRT7', 'CFTR', 'MUC1','BHLHA15', 'PTF1A', 'CPA1', 'AMY2A']:
  sc.pl.umap(adata, color= g, use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  plt.savefig("{0}/{1}_marker_gene_{2}_UMAPs.pdf".format(qcDir, bname,g) , bbox_inches='tight'); plt.close('all')

# Save the coordinates of the umap
rowIdx     = adata.obs.index.values.tolist()
umapCordDF = pd.DataFrame(adata.obsm['X_umap'],index=rowIdx, columns=['UMAP1','UMAP2'])
umapCordDF.to_csv("{0}/02_normalizedRaw_{1}_UMAP_Coordinated.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="PatId")

# Dataframe of ranked genes
for g in ['beta','alpha','ductal','acinar','delta','activated_stellate','gamma','endothelial','quiescent_stellate','macrophage','mast','epsilon','schwann','t_cell']:
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(adata.uns['rank_genes_Clusters'][n])[g]
  # Save dataframes
  ngDF.to_csv("{0}/02_normalizedRaw_{1}_rank_genes_Clusters_{2}.txt".format(countsDir, projName,g), sep='\t', header=True, index=False)

namesDF    = pd.DataFrame(adata.uns['rank_genes_Clusters']['names'])
namesDF.to_csv("{0}/02_normalizedRaw_{1}_rank_genes_Clusters_names.txt".format(countsDir, projName), sep='\t', header=True, index=False)

scoresDF   = pd.DataFrame(adata.uns['rank_genes_Clusters']['scores'])
scoresDF.to_csv("{0}/02_normalizedRaw_{1}_rank_genes_Clusters_scores.txt".format(countsDir, projName), sep='\t', header=True, index=False)

logfcDF    = pd.DataFrame(adata.uns['rank_genes_Clusters']['logfoldchanges'])
logfcDF.to_csv("{0}/02_normalizedRaw_{1}_rank_genes_Clusters_logfoldchanges.txt".format(countsDir, projName), sep='\t', header=True, index=False)

pvalsDF    = pd.DataFrame(adata.uns['rank_genes_Clusters']['pvals'])
pvalsDF.to_csv("{0}/02_normalizedRaw_{1}_rank_genes_Clusters_pvals.txt".format(countsDir, projName), sep='\t', header=True, index=False)

pvalsadjDF = pd.DataFrame(adata.uns['rank_genes_Clusters']['pvals_adj'])
pvalsadjDF.to_csv("{0}/02_normalizedRaw_{1}_rank_genes_Clusters_pvals_adj.txt".format(countsDir, projName), sep='\t', header=True, index=False)

