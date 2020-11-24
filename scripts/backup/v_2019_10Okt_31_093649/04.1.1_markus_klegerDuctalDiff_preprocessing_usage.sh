# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

# Transpose the raw genesymbols file for genes
datamash -H transpose < /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/input/scMatrix.txt > /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/input/T_scMatrix.txt
datamash -H transpose < /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/input/input_klegerDuctalDiff.txt > /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/input/T_input_klegerDuctalDiff.txt

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
projName        = "klegerDuctalDiff" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/input/01_preprocessing"; create_dir("{0}".format(output_dir))
# cc_genes_file   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/Macosko_cell_cycle_genes.txt"
minGenesPerCell = 100
minCellsPergene = 2
bname           = projName
qcDir           = "{0}/qc".format(output_dir); create_dir(qcDir)
countsDir       = "{0}/counts".format(output_dir); create_dir(countsDir)
# annotation_file = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/input/targets.txt"
input_matrix_file = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/klegerDuctalDiff/input/scMatrix.txt'

# Add 
# 1) Reading the data
orgdataDF   = pd.read_csv(input_matrix_file, sep="\t", index_col=['PatID'])
orgcountsDF = orgdataDF.copy()
orgcountsDF.drop(['DonorAge', 'Celltype'], axis=1, inplace=True)
adata = sc.AnnData(orgcountsDF)

# Make variable names unique
adata.var_names_make_unique()

# Add annotation
adata.obs['DonorAge'] = orgdataDF['DonorAge']
adata.obs['Celltype'] = orgdataDF['Celltype']

# # Convert the sparse count matrices to dense represntation
# adata.X = adata.X.toarray()
# AnnData object with n_obs × n_vars = 2527 × 23359 
#     obs: 'DonorAge', 'Celltype'


# Checking the total size of the data set
adata.shape # We have 2527 cells and 23359 genes in the dataset

# 2) Quality control 
# 2.1) Calculate QC covariates
adata.obs['n_counts'] = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes'] = (adata.X > 0).sum(1)

mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']

# 2.2) Plot QC metrics
# Sample quality plots
t1 = sc.pl.violin(adata, 'n_counts', groupby='DonorAge', size=2, log=True, cut=0, show=False)
plt.savefig("{0}/{1}_DonorAge_nCounts_plot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
t2 = sc.pl.violin(adata, 'mt_frac', groupby='DonorAge', show=False)
plt.savefig("{0}/{1}_DonorAge_mtFraction_plot.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

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

# Total number of cells: 2527
# Number of cells after min count filter: 2527
# Number of cells after MT filter: 2527
# Trying to set attribute `.obs` of view, making a copy.
# Number of cells after gene filter: 2527


# 2.6) Filter genes according to identified QC thresholds:
# Min 5 cells - filters out 0 count genes
print('Total number of genes: {:d}'.format(adata.n_vars))
sc.pp.filter_genes(adata, min_cells=10)
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
# Total number of genes: 23359
# Number of genes after cell filter: 18568

# 2.7) Save the filtered raw data as tab separated matrix 
adata.to_df().to_csv("{0}/01_raw_T_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="CellId")
adata.to_df().T.to_csv("{0}/01_raw_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol")

# 2.8) Calculations for the visualizations
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(adata, random_state = 2105)
sc.tl.umap(adata, random_state = 2105)

# 2.9) Plot visualizations
sc.pl.pca_scatter(adata, color='n_counts', show=False)
plt.savefig("{0}/01_raw_{1}_clustering_ncounts_PCA.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['DonorAge'], show=False)
plt.savefig("{0}/01_raw_{1}_clustering_DonorAge_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['Celltype'], show=False)
plt.savefig("{0}/01_raw_{1}_clustering_Celltype_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

########################
# 3) Expression recovery (denoising), Normalization and log transformation
# Denoise the data on raw counts
# sce.pp.dca(adata)
# NOTE: The denoising is causing the batches effects to get stronger and even with combat or BBKNm they are still are there with high impact.
#       For this Manec dataset, I`ll turn off the DCA
# # Save the denoised data
# adata.to_df().to_csv("{0}/02_denoised_DCA_T_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="CellId")
# adata.to_df().T.to_csv("{0}/02_denoised_DCA_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol")


# # Get the raw matrix to pass it along to Scnorm
# rawMat       = adata.to_df().T
# input_groups = adata.obs['DonorAge'].tolist()
# with open('conditions.txt', 'w') as f: f.write("\n".join(str(item) for item in input_groups))
# # Run Scnorm in R
# %%R -i rawMat -i input_groups -o normRawData
# groups <- (as.vector(unlist(input_groups)))
# normRawMat  <- SCnorm(Data=rawMat, Conditions=groups, K=2)
# normRawData <- results(normRawMat, type = "NormalizedData")
# # Keep the count data in a counts layer
# adata.layers["counts"] = adata.X.copy()
# # Normalize and log adata 
# adata.X = normRawData
# sc.pp.log1p(adata)

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
sc.pp.combat(adata, key='DonorAge')

# 4.2) Highly Variable Genes
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig("{0}/{1}_highly_variable_genes.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Highly variable gene information is stored automatically in the adata.var['highly_variable'] field. The dataset now contains:
# a 'counts' layer with count data
# log-normalized data in adata.raw
# batch corrected data in adata.X
# highly variable gene annotations in adata.var['highly_variable']
# The HVG labels will be used to subselect genes for clustering and trajectory analysis.

# # 5) Biological correction
# # 5.1) Cell cycle scoring
# #Score cell cycle and visualize the effect:
# cc_genes         = pd.read_table(cc_genes_file, delimiter='\t')
# s_genes          = cc_genes['S'].dropna()
# g2m_genes        = cc_genes['G2.M'].dropna()
# s_genes_mm       = [gene.lower().capitalize() for gene in s_genes]
# g2m_genes_mm     = [gene.lower().capitalize() for gene in g2m_genes]
# s_genes_mm_ens   = adata.var_names[np.in1d(adata.var_names, s_genes_mm)]
# g2m_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, g2m_genes_mm)]

# sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)
# sc.pl.umap(adata, color=['S_score', 'G2M_score'], use_raw=False, show=False)
# plt.savefig("{0}/{1}_S_G2M_Phase_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
# sc.pl.umap(adata, color='phase', use_raw=False, show=False)
# plt.savefig("{0}/{1}_combined_Phase_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# Save the normalized, log transformed, batch and cell cycle corrected data
adata.to_df().to_csv("{0}/02_normalizedRaw_T_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="PatId")
adata.to_df().T.to_csv("{0}/02_normalizedRaw_{1}_filtered.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol")

# 7) Clustering
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
adata.obs['louvain'].value_counts()
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
sc.pl.umap(adata, color=['louvain','louvain_r1', 'louvain_r0.5'], palette=sc.pl.palettes.default_64, show=False)
plt.savefig("{0}/{1}_clustering_louvain_r1_r05_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# UMAPs
sc.pl.pca_scatter(adata, color='n_counts', show=False)
plt.savefig("{0}/02_norm_{1}_clustering_ncounts_PCA.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['DonorAge'], show=False)
plt.savefig("{0}/02_norm_{1}_clustering_DonorAge_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['Celltype'], show=False)
plt.savefig("{0}/02_norm_{1}_clustering_Celltype_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color='n_counts', show=False)
plt.savefig("{0}/02_norm_{1}_clustering_ncounts_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')
sc.pl.umap(adata, color=['log_counts', 'mt_frac'], show=False)
plt.savefig("{0}/02_norm_{1}_clustering_logCounts_mtFrac_UMAP.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

# 7.2) Marker genes & cluster annotation
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby='louvain_r0.5', key_added='rank_genes_r0.5')
sc.tl.rank_genes_groups(adata, groupby='louvain_r1', key_added='rank_genes_r1')
sc.tl.rank_genes_groups(adata, groupby='louvain', key_added='rank_genes')

# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes', groups=['0','1','2'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_0_1_2.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.rank_genes_groups(adata, key='rank_genes', groups=['3','4','5'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_3_4_5.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.rank_genes_groups(adata, key='rank_genes', groups=['6', '7', '8'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_7_6_8.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.rank_genes_groups(adata, key='rank_genes', groups=['9', '10', '11'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_7_6_8.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.rank_genes_groups(adata, key='rank_genes', groups=['12', '13', '14'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_7_6_8.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')


# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

sc.pl.umap(adata, color=['louvain', 'KRT19'], use_raw=False, color_map=mymap)

















plt.savefig("{0}/{1}_marker_genes_adult_stomach_mouse_cell_atlas_UMAPs.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')


sc.pl.umap(adata, color=['Ctrb1','Cd2', 'Myc', 'Rb1','Cdkn2a'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)
sc.pl.umap(adata, color=['Birc5', 'Casp3', 'Stat3', 'Alb'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

sc.pl.umap(adata, color=['Birc5', 'Stat3', 'Reg1', 'Gm26917', 'Ctrb1', 'Clps', 'Hbb-bs', 'Hba-a1','Hba-a2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)
sc.pl.umap(adata, color=['Birc5', 'Blvrb', 'Car2', 'Hbb-bt', 'Clps', 'Hbb-bs', 'Hba-a1','Hba-a2', 'Hspa1b', 'Apoe', 'C1qb', 'Cxcl2', 'Slpi'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)
sc.pl.umap(adata, color=['Sh2d1a', 'Cd3d','Cd3e','Cd8a','Retnlg','S100a8','S100a9','Cxcl2', 'Slpi', 'Srgn', 'Cd84', 'Stip1','Cd44', 'Jak1'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

sc.pl.umap(adata, color=['Btg1', 'Morf4l2','Marcks','Ptprc','BC005537','Ctnnb1','Ptma','AW112010', 'Hnrnpf', 'Hspa1b', 'Hnrnph1', 'Dazap2','Laptm5', 'Id2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

sc.pl.umap(adata, color=['louvain_r0.5','Apoe', 'Dcn', 'Cd74','Pf4', 'Lyz2', 'Ly6c1', 'Plvap', 'Fabp5', 'Birc5', 'Ube2c', 'Dmbt1', 'Cyr61', 'Igfbp3', 'Clps'], use_raw=False, color_map='hot', show=False)
plt.savefig("{0}/{1}_marker_genes_adult_stomach_mouse_cell_atlas_UMAPs.png".format(qcDir, bname) , bbox_inches='tight'); plt.close('all')


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
  plt.savefig("{0}/{1}_marker_genes_stomach_{2}_violinPlots.png".format(qcDir, bname, k) , bbox_inches='tight'); plt.close('all')
  sc.pl.umap(adata, color=['louvain','{0}_marker_expr'.format(k)], color_map=mymap, show=False)
  plt.savefig("{0}/{1}_marker_genes_stomach_{2}_UMAPs.png".format(qcDir, bname, k) , bbox_inches='tight'); plt.close('all')
  

# Generate the UMAPs for each marker categories
plt.figure(figsize=(100,8))
for k in ma_marker_genes.keys():
  ids = np.in1d(adata.var_names, ma_marker_genes[k])
  adata.obs['{0}_ma_marker_expr'.format(k)] = adata.X[:,ids].mean(1)
  # sc.pl.violin(adata, '{0}_ma_marker_expr'.format(k), groupby='louvain_r0.5', show=False)
  # plt.savefig("{0}/{1}_mouse_cellatlas_marker_genes_stomach_{2}_violinPlots.png".format(qcDir, bname, k) , bbox_inches='tight'); plt.close('all')
  # sc.pl.umap(adata, color=['louvain','{0}_ma_marker_expr'.format(k)], color_map=mymap, show=False)
  sc.pl.umap(adata, color=['{0}_ma_marker_expr'.format(k)], color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  plt.savefig("{0}/{1}_mouse_cellatlas_marker_genes_stomach_{2}_UMAPs.png".format(qcDir, bname, k) , bbox_inches='tight', dpi=300); plt.close('all')
  
# Categories to rename
adata.obs['louvain_r0.5'].cat.categories

# On the subclusters
sbadata =adata.copy()
# sc.pl.umap(sbadata, color='louvain', palette=sc.pl.palettes.vega_20, size=50)
plt.figure(figsize=(100,8))
sc.pl.umap(sbadata, color='louvain_r1', palette=sc.pl.palettes.vega_20, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/{1}_Louvainr05_UMAPs.png".format(qcDir, bname, k) , bbox_inches='tight', dpi=300); plt.close('all')

# Groups interested = 0, 2, 5, 10
sc.tl.rank_genes_groups(sbadata, groupby='louvain_r1', key_added='rank_genes_0_2', groups = '2', reference = '0')
sc.pl.rank_genes_groups(sbadata, key='rank_genes_0_2', fontsize=12)


grp = ['4', '6', '13']
ref = '11'
sc.tl.rank_genes_groups(adata, groupby='louvain_r1', key_added='rank_genes_{0}_{1}'.format(grp, ref), groups = list(grp), reference =  '{}'.format(ref))
sc.pl.rank_genes_groups(adata, key='rank_genes_{0}_{1}'.format(grp, ref), fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_pairwise_{2}_{3}.png".format(qcDir, bname, grp, ref) , bbox_inches='tight'); plt.close('all')

# Plot umap for the top genes for each cluster
clusterDir = "{0}/clusterDir".format(qcDir); create_dir(clusterDir)
for cid in [0, 2, 5, 10]:
  cluster = adata[(adata.obs['louvain_r1'] == '{0}'.format(cid))].to_df() 
  topGenes = np.mean(cluster, axis=0).sort_values(ascending=False)[0:25].index.tolist()
  print("\n- Cluster{0}: {1}".format(cid, topGenes))
  plt.figure(figsize=(100,8))
  sc.pl.umap(adata, color=topGenes, use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  plt.savefig("{0}/{1}_top_25_genes_cluster{2}_UMAPs.png".format(clusterDir, bname, cid) , bbox_inches='tight'); plt.close('all')
  
# Plot umap for the top marker genes for each cluster
markerDir = "{0}/markerDir/pdf".format(qcDir); create_dir(markerDir)
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
plt.savefig("{0}/{1}_Clusters_CellTypes_UMAPs.png".format(qcDir, bname, k) , bbox_inches='tight', dpi=300); plt.close('all')


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