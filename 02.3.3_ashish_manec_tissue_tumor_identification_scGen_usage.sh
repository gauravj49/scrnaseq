# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

# Source: https://nbviewer.jupyter.org/github/M0hammadL/scGen_notebooks/blob/master/notebooks/scgen_batch_removal.ipynb

ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)

# Loading the python libraries
import scanpy as sc
import scanpy.external as sce
from gjainPyLib import *
import pickle
import scgen

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

# System variables and directories
projName        = "manec_tissues_merged_except1079" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/tissue/02_scGenMitTumorId/{0}".format(projName); create_dir("{0}".format(output_dir))
cc_genes_file   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/Macosko_cell_cycle_genes.txt"
minGenesPerCell = 100
minCountPerCell = 300
minCellsPergene = 25
bname           = projName
plotsDir        = "{0}/plots".format(output_dir); create_dir(plotsDir)
dataDir         = "{0}/data".format(output_dir) ; create_dir(dataDir)

# 1 Reading the data
# Merge 10x datasets for different mices
# https://github.com/theislab/scanpy/issues/267
#                     'input/manec/pilot2/bulk1079_mouse_filtered_feature_bc_matrix.h5'
tissueFilenames = [
                    'input/manec/tissue/bulk997_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/bulk1001_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/bulk1018_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/stomach1001_mouse_filtered_feature_bc_matrix.h5'
                  ]
adatas          = [sc.read_10x_h5(f) for f in tissueFilenames]
adata           = adatas[0].concatenate(adatas[1:])

# Make variable names unique
adata.var_names_make_unique()

# Add tissue id column for the batches
tissueIdDict =  {
                    '0':'bulk997', 
                    '1':'bulk1001', 
                    '2':'bulk1018', 
                    '3':'stomach1001'
                }
adata.obs['tissueID'] = adata.obs['batch'].map(tissueIdDict)

# Convert the sparse count matrices to dense represntation
adata.X = adata.X.toarray()

# Checking the total size of the data set
adata.shape # We have 7987 cells and 55488 genes in the dataset

# Make a copy of the original data
origadata = adata.copy()

# 2) Quality control 
# 2.1) Calculate QC covariates
adata.obs['n_counts']   = adata.X.sum(1)
adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
adata.obs['n_genes']    = (adata.X > 0).sum(1)
mt_gene_mask            = [gene.startswith('mt-') for gene in adata.var_names]
adata.obs['mt_frac']    = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']

# 2.2) Plot QC metrics
# Sample quality plots
t1 = sc.pl.violin(adata, 'n_counts', groupby='tissueID', size=2, log=True, cut=0, show=False)
plt.savefig("{0}/01_{1}_tissueID_nCounts_plot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
t2 = sc.pl.violin(adata, 'mt_frac', groupby='tissueID', show=False)
plt.savefig("{0}/01_{1}_tissueID_mtFraction_plot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.3) Data quality summary plots
p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/01_{1}_genes_counts_scatterplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

p2 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/01_{1}_genes_counts_scatterplot_zoomedin.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.3) Thresholding decision based on counts
p3 = sns.distplot(adata.obs['n_counts'], kde=False); #plt.show()
plt.savefig("{0}/01_{1}_ncounts_histogramplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p4 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<2000], kde=False, bins=1000); #plt.show()
plt.savefig("{0}/01_{1}_ncounts_histogramplot_lessthan_2000.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p5 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>5000], kde=False, bins=1000); #plt.show()
plt.savefig("{0}/01_{1}_ncounts_histogramplot_greaterthan_5000.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.4) Thresholding decision based on genes
p6 = sns.distplot(adata.obs['n_genes'], kde=False, bins=1000); # plt.show()
plt.savefig("{0}/01_{1}_genes_histogramplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p7 = sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, bins=500); # plt.show()
plt.savefig("{0}/01_{1}_genes_histogramplot_lessthan_1000.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.5) Filter cells according to identified QC thresholds:
origadata = adata.copy()
print('Total number of cells: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_counts = minCountPerCell)
print('Number of cells after min count filter: {:d}'.format(adata.n_obs))

# sc.pp.filter_cells(adata, max_counts = 40000)
# print('Number of cells after max count filter: {:d}'.format(adata.n_obs))

adata = adata[adata.obs['mt_frac'] < 0.25]
print('Number of cells after MT filter: {:d}'.format(adata.n_obs))

sc.pp.filter_cells(adata, min_genes = minGenesPerCell)
print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

# Total number of cells: 7987
# Number of cells after min count filter: 7987
# Number of cells after MT filter: 7577
# Trying to set attribute `.obs` of view, making a copy.
# Number of cells after gene filter: 7577

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
sc.pl.umap(adata, color=['tissueID'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/01_raw_{1}_clustering_tissueID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adata, color=['tissueID'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/01_raw_{1}_clustering_tissueID_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

########################
# 3) Normalization using scran
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
plt.savefig("{0}/02_{1}_scrna_sizefactors_vs_ncounts.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.scatter(adata, 'size_factors', 'n_genes' , show=False)
plt.savefig("{0}/02_{1}_scrna_sizefactors_vs_ngenes.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sns.distplot(size_factors, bins=50, kde=False)
plt.savefig("{0}/02_{1}_scrna_sizefactors_histogram.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Keep the count data in a counts layer
adata.layers["counts"] = adata.X.copy()

# 3.5) Normalize adata 
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)

# Store the full data set in 'raw' as log-normalised data for statistical testing
adata.raw = adata

#######################
# Add celltypes
cellTypesDF = pd.read_csv('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/tissue/01_mitTumorIdentification/manec_tissues_merged_except1079/counts/03_manec_tissues_merged_except1079_cellTypes.txt', sep="\t", index_col=0)
adata.obs['cellType'] = cellTypesDF['cellType']
########################
# 4) Batch correction using scGen
# 4.1) Creating and training the model
network = scgen.VAEArith(x_dimension= adata.shape[1], model_path="./models/batch" )

# 4.2) Train the model for 100 epochs
network.train(train_data=adata, n_epochs=500)

# 4.3) Batch Removal
# Now, we pass our train (uncorected) adata and network object to batch_removal fucntion which returns adata of batch corrected data (cadata)
cadata =  scgen.batch_removal(network, adata, batch_key="tissueID", cell_label_key="cellType")

# Plot the UMAPs
sc.pp.neighbors(cadata)
sc.tl.umap(cadata, random_state = 2105, n_components=3)
# sc.pl.umap(cadata, color=["tissueID", "cellType"], wspace=.5)

sc.pl.umap(cadata, color=['tissueID'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_batchCorrected_{1}_clustering_tissueID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(cadata, color=['tissueID'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False, projection='3d')
plt.savefig("{0}/03_batchCorrected_{1}_clustering_tissueID_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.umap(cadata, color=['cellType'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_batchCorrected_{1}_clustering_cellType_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(cadata, color=['cellType'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False, projection='3d')
plt.savefig("{0}/03_batchCorrected_{1}_clustering_cellType_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.umap(cadata, color=['cellType'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False,legend_loc='on data')
plt.savefig("{0}/03_batchCorrected_{1}_clustering_cellType_legend_onData_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')


#########################################
# Save session
import dill
filename = "{0}/{1}.pkl".format(output_dir, projName)

# and to load the session again:
import dill
filename = "{0}/{1}.pkl".format(output_dir, projName)
dill.load_session(filename)
#########################################
