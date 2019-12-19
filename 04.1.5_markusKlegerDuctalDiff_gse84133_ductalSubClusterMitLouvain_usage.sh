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


















# Process subcluster 8 individually (but not preprocessing from the beginning)
# Get the cluster 8 anndata
adataCluster8 = adata[adata.obs['louvain_r1']=='8',:]

# Mucin and CFTR enriched ductal subpopulation
sc.pl.umap(adataCluster8, color=['MUC1', 'CFTR', 'TFF1', 'CD44'], use_raw=False, color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, ncols=2, legend_loc='on data', show=False)
plt.savefig("{0}/04_{1}_marker_genes_mucin_cftr_enriched_ductal_subpopulation_UMAPs.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adataCluster8, color=['MUC1', 'CFTR', 'TFF1', 'CD44'], use_raw=False, color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, ncols=2, legend_loc='on data', projection='3d', show=False)
plt.savefig("{0}/04_{1}_marker_genes_mucin_cftr_enriched_ductal_subpopulation_3D_UMAPs.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# Analyze ductal (1077 cells) subgroups (CFTR_high/MUC1_low and CFTR_low/MUC1_high) 
# Get the index of ductal subcluster
Cluster8ductalDF = adataCluster8[adataCluster8.obs['Cluster']=='ductal'].to_df()   # 1077

# Get the ratio of expression
# pd.set_option('display.float_format', lambda x: '%.2f' % x)
Cluster8ductalDF['CFTR_rank'] = Cluster8ductalDF['CFTR'].rank()
Cluster8ductalDF['MUC1_rank'] = Cluster8ductalDF['MUC1'].rank()
Cluster8ductalDF['rank_CFTR_minus_MUC1'] = Cluster8ductalDF.apply(lambda row: (row['CFTR_rank'] - row['MUC1_rank']), axis=1)
# Cluster8ductalDF[['CFTR', 'MUC1', 'rank_CFTR_minus_MUC1']]

# # Save the relevant data into an excel
# Cluster8ductalDF[['CFTR', 'MUC1', 'CFTR_rank', 'MUC1_rank', 'rank_CFTR_minus_MUC1']].to_csv("{0}/04_{1}_clean_ductal_subcluster_cluster8_gene_expression_ranks.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="CellId", float_format='%.4g')

# Get median for both genes
# cftrMedianExp = np.float(Cluster8ductalDF['CFTR'].describe()[['50%']].values)
# muc1MedianExp = np.float(Cluster8ductalDF['MUC1'].describe()[['50%']].values)
chmlDuctalIdx  = Cluster8ductalDF[Cluster8ductalDF['rank_CFTR_minus_MUC1'] >= 100 ].index.tolist()                # 105
clmhDuctalIdx  = Cluster8ductalDF[Cluster8ductalDF['rank_CFTR_minus_MUC1'] <= -100].index.tolist()                # 94
otherDuctalIdx = Cluster8ductalDF.loc[~Cluster8ductalDF.index.isin(chmlDuctalIdx + clmhDuctalIdx)].index.tolist() # 237

# print(len(chmlDuctalIdx))
# print(len(clmhDuctalIdx))
# print(len(otherDuctalIdx))

# Get a new subcluster column
adataCluster8.obs['ductal_subcluster_cluster8'] = adataCluster8.obs['Cluster']

# Add new categories
adataCluster8.obs['ductal_subcluster_cluster8'].cat.add_categories(['ductal_cftrHigh_muc1Low','ductal_cftrLow_muc1High', 'ductal_other'], inplace=True) 

# Get a new subcluster column
adataCluster8.obs['ductal_subcluster_cluster8'].loc[chmlDuctalIdx]  = adataCluster8.obs['ductal_subcluster_cluster8'].loc[chmlDuctalIdx].replace('ductal','ductal_cftrHigh_muc1Low') 
adataCluster8.obs['ductal_subcluster_cluster8'].loc[clmhDuctalIdx]  = adataCluster8.obs['ductal_subcluster_cluster8'].loc[clmhDuctalIdx].replace('ductal','ductal_cftrLow_muc1High') 
adataCluster8.obs['ductal_subcluster_cluster8'].loc[otherDuctalIdx] = adataCluster8.obs['ductal_subcluster_cluster8'].loc[otherDuctalIdx].replace('ductal','ductal_other') 

adataCleanDucCluster8 = adataCluster8[~((adataCluster8.obs['ductal_subcluster_cluster8'] == 'schwann') | (adataCluster8.obs['ductal_subcluster_cluster8'] == 'acinar') | (adataCluster8.obs['ductal_subcluster_cluster8'] == 'activated_stellate') | (adataCluster8.obs['ductal_subcluster_cluster8'] == 'ductal_other'))]

# Visualize new subclusers
sc.pl.umap(adataCluster8, color=['ductal_subcluster_cluster8'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/04_norm_{1}_ductal_subcluster_cluster8_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adataCluster8, color=['ductal_subcluster_cluster8'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/04_norm_{1}_ductal_subcluster_cluster8_UMAP_3D.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

sc.pl.umap(adataCleanDucCluster8, color=['ductal_subcluster_cluster8'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/04_norm_{1}_clean_ductal_subcluster_cluster8_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adataCleanDucCluster8, color=['ductal_subcluster_cluster8'], palette=sc.pl.palettes.vega_10, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/04_norm_{1}__clean_ductal_subcluster_cluster8_UMAP_UMAP_3D.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# Mucin and CFTR enriched in clean ductal subpopulation
sc.pl.umap(adataCleanDucCluster8, color=['MUC1', 'CFTR', 'TFF1', 'CD44'], use_raw=False, color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, ncols=2, legend_loc='on data', show=False)
plt.savefig("{0}/04_{1}_marker_genes_mucin_cftr_enriched_clean_ductal_subpopulation_UMAPs.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adataCleanDucCluster8, color=['MUC1', 'CFTR', 'TFF1', 'CD44'], use_raw=False, color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, ncols=2, legend_loc='on data', projection='3d', show=False)
plt.savefig("{0}/04_{1}_marker_genes_mucin_cftr_enriched_clean_ductal_subpopulation_3D_UMAPs.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')


# 9.1) Marker genes & cluster annotation
# Calculate marker genes
sc.tl.rank_genes_groups(adataCleanDucCluster8, groupby='ductal_subcluster_cluster8', groups=['ductal_cftrHigh_muc1Low'], key_added='rank_genes_CftrHighMuc1Low_over_Muc1HighCftrLow', reference='ductal_cftrLow_muc1High', n_genes=adataCleanDucCluster8.shape[1])
plt.savefig("{0}/04_{1}_marker_genes_ranking_rank_genes_CftrHighMuc1Low_over_Muc1HighCftrLow.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Extract and save dataframe of ranked genes
ngDF = pd.DataFrame()
for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
  ngDF[n] = pd.DataFrame(adataCleanDucCluster8.uns['rank_genes_CftrHighMuc1Low_over_Muc1HighCftrLow'][n])['ductal_cftrHigh_muc1Low']
ngDF.to_csv("{0}/03_normalizedRaw_CleanDuctalSubCluster_rank_genes_CftrHighMuc1Low_over_Muc1HighCftrLow.txt".format(countsDir), sep='\t', header=True, index=False, float_format='%.2g')

# Plot additional marker genes on UMAP
# Mucin enriched ductal subpopulation
sc.pl.umap(adataCleanDucCluster8, color=['ductal_subcluster_cluster8', 'MUC1', 'MUC13', 'TFF1', 'TFF2'], use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/03_{1}_marker_genes_mucin_enriched_ductal_subpopulation_UMAPs_3D.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# CFTR population
sc.pl.umap(adataCleanDucCluster8, color=['ductal_subcluster_cluster8', 'KRT19', 'KRT7', 'CFTR', 'AQP3', 'AQP5', 'CA2', 'CA4', 'SCTR', 'SLC26A6'], use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/03_{1}_marker_genes_cftr_enriched_ductal_subpopulation_UMAPs_3D.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

# CFTR population
sc.pl.umap(adataCleanDucCluster8, color=['ductal_subcluster_cluster8', 'ANXA2','ANXA3','ANXA4','CDX2','ID2','SOX6','SCTR','CD44','ID3','CLDN4','KLF6','LGALS3','NCOA7','ANXA1','ANXA2','CLDN1','KRT17','LAMC2','S100A16','ABCC3','AMBP','ANXA4','CLDN18','ONECUT2','PDX1','NKX6.1'], use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/03_{1}_additional_interesting_genes_UMAPs_3D.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')

for g in ['ductal_subcluster_cluster8', 'MUC1', 'MUC13', 'TFF1', 'TFF2', 'KRT19', 'KRT7', 'CFTR', 'AQP3', 'AQP5', 'CA2', 'CA4', 'SCTR', 'SLC26A6', 'ANXA2','ANXA3','ANXA4','CDX2','ID2','SOX6','SCTR','CD44','ID3','CLDN4','KLF6','LGALS3','NCOA7','ANXA1','ANXA2','CLDN1','KRT17','LAMC2','S100A16','ABCC3','AMBP','ANXA4','CLDN18','ONECUT2','PDX1','NKX6.1']:
  sc.pl.umap(adataCleanDucCluster8, color= g, use_raw=False, color_map=mymap, size=50, edgecolor='grey', linewidth=0.01, alpha=0.9, projection='3d', show=False)
  plt.savefig("{0}/03_{1}_marker_gene_{2}_UMAPs_3D.pdf".format(qcDir, bname,g) , bbox_inches='tight'); plt.close('all')


###############################################################################################################################
# Analysis of Cluster 7 + Cluster 8

# Get the cluster 7+8 anndata
adataCluster78 = adata[(adata.obs['louvain_r1']=='7')|(adata.obs['louvain_r1']=='8'),:]

# Mucin and CFTR enriched ductal subpopulation
sc.pl.umap(adataCluster78, color=['MUC1', 'CFTR', 'TFF1', 'CD44'], use_raw=False, color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, ncols=2, legend_loc='on data', show=False)
plt.savefig("{0}/04_{1}_marker_genes_mucin_cftr_enriched_ductal_subpopulation_UMAPs.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')
sc.pl.umap(adataCluster78, color=['MUC1', 'CFTR', 'TFF1', 'CD44'], use_raw=False, color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, ncols=2, legend_loc='on data', projection='3d', show=False)
plt.savefig("{0}/04_{1}_marker_genes_mucin_cftr_enriched_ductal_subpopulation_3D_UMAPs.png".format(qcDir, bname) , bbox_inches='tight', dpi=300); plt.close('all')


# Analyze ductal (1077 cells) subgroups (CFTR_high/MUC1_low and CFTR_low/MUC1_high) 
# Get the index of ductal subcluster
Cluster8ductalDF = adataCluster8[adataCluster8.obs['Cluster']=='ductal'].to_df()   # 1077

# Get the ratio of expression
# pd.set_option('display.float_format', lambda x: '%.2f' % x)
Cluster8ductalDF['CFTR_rank'] = Cluster8ductalDF['CFTR'].rank()
Cluster8ductalDF['MUC1_rank'] = Cluster8ductalDF['MUC1'].rank()
Cluster8ductalDF['rank_CFTR_minus_MUC1'] = Cluster8ductalDF.apply(lambda row: (row['CFTR_rank'] - row['MUC1_rank']), axis=1)
# Cluster8ductalDF[['CFTR', 'MUC1', 'rank_CFTR_minus_MUC1']]

# Get median for both genes
# cftrMedianExp = np.float(Cluster8ductalDF['CFTR'].describe()[['50%']].values)
# muc1MedianExp = np.float(Cluster8ductalDF['MUC1'].describe()[['50%']].values)
chmlDuctalIdx  = Cluster8ductalDF[Cluster8ductalDF['rank_CFTR_minus_MUC1'] >= 100 ].index.tolist()                # 105
clmhDuctalIdx  = Cluster8ductalDF[Cluster8ductalDF['rank_CFTR_minus_MUC1'] <= -100].index.tolist()                # 94
otherDuctalIdx = Cluster8ductalDF.loc[~Cluster8ductalDF.index.isin(chmlDuctalIdx + clmhDuctalIdx)].index.tolist() # 237

# print(len(chmlDuctalIdx))
# print(len(clmhDuctalIdx))
# print(len(otherDuctalIdx))

# Get a new subcluster column
adataCluster8.obs['ductal_subcluster_cluster8'] = adataCluster8.obs['Cluster']

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
