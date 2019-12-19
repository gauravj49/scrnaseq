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
projName          = "klegerDuctalDiff_gse84133" # MANEC_merged_except1079_hMYC_forcecells
input_matrix_file = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/02_gse84133/scMatrix.txt'
output_dir        = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/02_gse84133"; create_dir("{0}".format(output_dir))
minGenesPerCell   = 100
minCellsPergene   = 2
bname             = projName
qcDir             = "{0}/analysis".format(output_dir); create_dir(qcDir)
countsDir         = "{0}/counts".format(output_dir); create_dir(countsDir)

# Define a nice colour map for gene expression
colors2    = plt.cm.Reds(np.linspace(0, 1, 128))
colors3    = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap      = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

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
origadata   = sc.AnnData(orgcountsDF)

# Work on a copy of original data
adata = origadata.copy()

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
plt.savefig("{0}/01_{1}_Batch_nCounts_plot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# 2.3) Thresholding decision based on counts
p3 = sns.distplot(adata.obs['n_counts'], kde=False); #plt.show()
plt.savefig("{0}/01_{1}_ncounts_histogramplot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
p4 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<2000], kde=False, bins=1000); #plt.show()
plt.savefig("{0}/01_{1}_ncounts_histogramplot_lessthan_2000.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
p5 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>5000], kde=False, bins=1000); #plt.show()
plt.savefig("{0}/01_{1}_ncounts_histogramplot_greaterthan_5000.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# 2.4) Thresholding decision based on genes
p6 = sns.distplot(adata.obs['n_genes'], kde=False, bins=1000); # plt.show()
plt.savefig("{0}/01_{1}_genes_histogramplot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# 2.5) Filter cells according to identified QC thresholds:
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
# adata.to_df().to_csv("{0}/01_raw_T_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="CellId", float_format='%.2g')
# adata.to_df().T.to_csv("{0}/01_raw_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol", float_format='%.2g')

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
plt.savefig("{0}/01_{1}_scrna_sizefactors_vs_ncounts.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.scatter(adata, 'size_factors', 'n_genes' , show=False)
plt.savefig("{0}/01_{1}_scrna_sizefactors_vs_ngenes.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sns.distplot(size_factors, bins=50, kde=False)
plt.savefig("{0}/01_{1}_scrna_sizefactors_histogram.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

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
plt.savefig("{0}/01_{1}_highly_variable_genes.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Save the normalized, log transformed, batch and cell cycle corrected data
normCorrectedDF            = adata.to_df()
normCorrectedDF['Batch']   = orgdataDF['Batch'].astype(str)
normCorrectedDF['Cluster'] = orgdataDF['Cluster']
# normCorrectedDF.to_csv("{0}/02_normalizedRaw_T_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="PatId", float_format='%.2g')
# normCorrectedDF.T.to_csv("{0}/02_normalizedRaw_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol", float_format='%.2g')

# 4.3) Visualizations
# Calculations for the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata, random_state = 2105, n_components=3)

# UMAPs for QC
sc.pl.pca_scatter(adata, color='n_counts', show=False)
plt.savefig("{0}/02_norm_{1}_ncounts_PCA.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.umap(adata, color=['Batch'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_Batch_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

sc.pl.umap(adata, color='n_counts', palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_ncounts_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# UMAPs with cluster information
# 2D plots
sc.pl.umap(adata, color=['Cluster'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_Cluster_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adata, color=['Cluster'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, legend_loc='on data', show=False)
plt.savefig("{0}/02_norm_{1}_Cluster_legendOnData_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# 3D plots
sc.pl.umap(adata, color=['Cluster'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/02_norm_{1}_Cluster_UMAP_3D.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# 5) LOUVAIN CLUSTERING 
# 5.1) Perform clustering - using highly variable genes
sc.tl.louvain(adata, resolution=1.0, key_added='louvain_r1'  , random_state=2105)
sc.pl.umap(adata, color=['louvain_r1'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, legend_loc='on data', show=False)
plt.savefig("{0}/02_norm_{1}_louvain_r1_legendOnData_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adata, color=['louvain_r1'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/02_norm_{1}_louvain_r1_UMAP_3D.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# Get a new subcluster column
adata.obs['customClusters'] = adata.obs['louvain_r1']

# Add new categories
adata.obs['customClusters'].cat.add_categories(['alpha','acinar','beta','delta', 'gamma','ductal'], inplace=True) 

# Get a new subcluster column
# 0+1+15   = alpha
# 2+10     = acinar
# 3+4+5+16 = beta
# 6        = delta
# 11       = gamma
# 7+8      = ductal
adata.obs['customClusters'].loc[adata.obs['customClusters']=='6']  = 'delta'
adata.obs['customClusters'].loc[adata.obs['customClusters']=='11'] = 'gamma'
adata.obs['customClusters'].loc[(adata.obs['customClusters']=='7')|(adata.obs['customClusters']=='8')]  = 'ductal'
adata.obs['customClusters'].loc[(adata.obs['customClusters']=='2')|(adata.obs['customClusters']=='10')] = 'acinar'
adata.obs['customClusters'].loc[(adata.obs['customClusters']=='0')|(adata.obs['customClusters']=='1')|(adata.obs['customClusters']=='15')] = 'alpha'
adata.obs['customClusters'].loc[(adata.obs['customClusters']=='3')|(adata.obs['customClusters']=='4')|(adata.obs['customClusters']=='5')|(adata.obs['customClusters']=='16')] = 'beta'

# Remove the cells that are not needed
# 9     = stellate
# 12    = endothelial
# 13,14 = unknown, might be alpha cells
# 17    = immune cells
# 18    = unknown cells
adataCustom = adata[~((adata.obs['customClusters']=='9') | (adata.obs['customClusters']=='12') |(adata.obs['customClusters']=='13') | (adata.obs['customClusters']=='14') | (adata.obs['customClusters']=='17') | (adata.obs['customClusters']=='18'))]

# Calculations for the visualizations
sc.pp.pca(adataCustom, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adataCustom)
sc.tl.umap(adataCustom, random_state = 2105, n_components=3)

# Mucin and CFTR enriched ductal subpopulation
sc.pl.umap(adataCustom, color=['MUC1', 'CFTR', 'TFF1', 'CD44'], use_raw=False, color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, ncols=2, legend_loc='on data', show=False)
plt.savefig("{0}/04_{1}_marker_genes_mucin_cftr_enriched_custom_subpopulation_UMAPs.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adataCustom, color=['MUC1', 'CFTR', 'TFF1', 'CD44'], use_raw=False, color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, ncols=2, legend_loc='on data', projection='3d', show=False)
plt.savefig("{0}/04_{1}_marker_genes_mucin_cftr_enriched_custom_subpopulation_3D_UMAPs.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# UMAPs with cluster information
# 2D plots
sc.pl.umap(adataCustom, color=['customClusters'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/04_{1}_CustomClusters_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adataCustom, color=['customClusters'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, legend_loc='on data', show=False)
plt.savefig("{0}/04_{1}_CustomClusters_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# 3D plots
sc.pl.umap(adataCustom, color=['customClusters'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/04_{1}_CustomClusters_UMAP_3D.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# 9.1) Marker genes analysis
sc.tl.rank_genes_groups(adataCustom, groupby='customClusters', key_added='rank_genes_customClusters', n_genes=adata.shape[1])

# adataCustom.obs['customClusters'].value_counts()                                                                                                                           Out[67]: 
# beta      2362
# alpha     2100
# acinar    1126
# ductal     910
# delta      607
# gamma      267
# Name: customClusters, dtype: int64

# Plot marker genes
sc.pl.rank_genes_groups(adataCustom, key='rank_genes_customClusters', groups=['beta','alpha','ductal','acinar','delta','gamma'], fontsize=12, show=False)
plt.savefig("{0}/04_{1}_marker_genes_ranking_customClusters.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Save UMAP of marker genes
sc.pl.umap(adataCustom, color=['customClusters', 'KRT19', 'KRT7', 'CFTR', 'MUC1','PTF1A', 'CPA1', 'AMY2A', 'CTRC', 'PRSS1'], use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/04_{1}_selected_marker_genes_customClusters_UMAPs.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
for g in ['customClusters', 'KRT19', 'KRT7', 'CFTR', 'MUC1', 'PTF1A', 'CPA1', 'AMY2A']:
  sc.pl.umap(adataCustom, color= g, use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  plt.savefig("{0}/04_{1}_marker_gene_{2}_UMAPs.pdf".format(qcDir, bname,g) , bbox_inches='tight'); plt.close('all')

# Save the coordinates of the umap
rowIdx     = adataCustom.obs.index.values.tolist()
umapCordDF = pd.DataFrame(adataCustom.obsm['X_umap'],index=rowIdx, columns=['UMAP1','UMAP2', 'UMAP3'])
umapCordDF.to_csv("{0}/04_normalizedRaw_{1}_customClusters_UMAP_Coordinates.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="PatId")

# Dataframe of ranked genes
for g in ['beta','alpha','ductal','acinar','delta','gamma']:
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(adataCustom.uns['rank_genes_customClusters'][n])[g]
  # Save dataframes
  ngDF.to_csv("{0}/04_normalizedRaw_{1}_rank_genes_customClusters_{2}.txt".format(countsDir, projName,g), sep='\t', header=True, index=False, float_format='%.2g')

namesDF    = pd.DataFrame(adataCustom.uns['rank_genes_customClusters']['names'])
namesDF.to_csv("{0}/04_normalizedRaw_{1}_rank_genes_customClusters_names.txt".format(countsDir, projName), sep='\t', header=True, index=False, float_format='%.2g')

scoresDF   = pd.DataFrame(adataCustom.uns['rank_genes_customClusters']['scores'])
scoresDF.to_csv("{0}/04_normalizedRaw_{1}_rank_genes_customClusters_scores.txt".format(countsDir, projName), sep='\t', header=True, index=False, float_format='%.2g')

logfcDF    = pd.DataFrame(adataCustom.uns['rank_genes_customClusters']['logfoldchanges'])
logfcDF.to_csv("{0}/04_normalizedRaw_{1}_rank_genes_customClusters_logfoldchanges.txt".format(countsDir, projName), sep='\t', header=True, index=False, float_format='%.2g')

pvalsDF    = pd.DataFrame(adataCustom.uns['rank_genes_customClusters']['pvals'])
pvalsDF.to_csv("{0}/04_normalizedRaw_{1}_rank_genes_customClusters_pvals.txt".format(countsDir, projName), sep='\t', header=True, index=False, float_format='%.2g')

pvalsadjDF = pd.DataFrame(adataCustom.uns['rank_genes_customClusters']['pvals_adj'])
pvalsadjDF.to_csv("{0}/04_normalizedRaw_{1}_rank_genes_customClusters_pvals_adj.txt".format(countsDir, projName), sep='\t', header=True, index=False, float_format='%.2g')

# Save the normalized, log transformed, batch and cell cycle corrected data
customClustersDF                              = adataCustom.to_df()
customClustersDF['originalClusterAssignment'] = adataCustom.obs['Cluster']
customClustersDF['Batch']                     = adataCustom.obs['Batch']
customClustersDF['louvain_r1']                = adataCustom.obs['louvain_r1']
customClustersDF['customClusters']            = adataCustom.obs['customClusters']
customClustersDF.to_csv("{0}/04_normalizedRaw_T_{1}_customClusters.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="PatId", float_format='%.2g')
customClustersDF.T.to_csv("{0}/04_normalizedRaw_{1}_customClusters.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol", float_format='%.2g')









#########################################
# Save session
import dill
filename = "{0}/{1}_0415_ductalSubClusterMitLouvain_usage.pkl".format(output_dir, projName)
dill.dump_session(filename)

# and to load the session again:
import dill
filename = "{0}/{1}_0415_ductalSubClusterMitLouvain_usage.pkl".format(output_dir, projName)
dill.load_session(filename)

#########################################
