# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

# Source: 
#   - https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb
#   - https://nbviewer.jupyter.org/github/theislab/trVAE/blob/master/examples/trVAE_Haber.ipynb


# NOTE: plt.subplot(number_of_rows, number_of_columns, index_of_the_subplot) 
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
import matplotlib.font_manager               # Incase of missing fonts (regenerate fonts)
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
projName        = "trvae_alltissues_except1079" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/{0}".format(projName); create_dir("{0}".format(output_dir))
ccGenes_macosko = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/macosko_cell_cycle_genes_mmu.txt"
ccGenes_regev   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/regev_lab_cell_cycle_genes_mmu.txt"
minGenesPerCell = 150
minCountPerCell = 100
maxCountPerCell = 50000 
minCellsPergene = 75
mtGenesFilter   = 0.25
rbGenesFilter   = 0.15
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
# Memory usage: current 1.02 GB, difference +1.02 GB
# Running Scanpy 1.4.5.1, on 2020-04-17 10:28.
# scanpy==1.4.5.1 anndata==0.7.1 umap==0.3.10 numpy==1.17.3 scipy==1.4.1 pandas==1.0.3 scikit-learn==0.21.3 statsmodels==0.10.1 python-igraph==0.7.1 louvain==0.6.1

# 1) Reading and performing QC on individual datasets

# 1.1) Reading the data in the anndata object individually
# Outlier_left_out: 'input/manec/pilot2/bulk1079_mouse_filtered_feature_bc_matrix.h5'
tissueFilenames = [
                    'input/manec/tissue/bulk997_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/bulk1001_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/bulk1018_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/stomach1001_mouse_filtered_feature_bc_matrix.h5'
                  ]
adatas          = [sc.read_10x_h5(f) for f in tissueFilenames]
adatas
# [AnnData object with n_obs × n_vars = 4334 × 55488 
#  AnnData object with n_obs × n_vars = 1852 × 55488 
#  AnnData object with n_obs × n_vars =  897 × 55488 
#  AnnData object with n_obs × n_vars =  904 × 55488]
#  var: 'gene_ids', 'feature_types', 'genome'

# Get the dictionary of tissue ids
tissueIdDict =  {
                    '0':'bulk997', 
                    '1':'bulk1001', 
                    '2':'bulk1018', 
                    '3':'stomach1001'
                }

# Create QC filtered adatas list
fadatas = list()

# 1.2) Make the variable names unique for each annData object separately and calculate some general qc-stats for genes and cells
minGenesPerCell = 5
minCountPerCell = 50
maxCountPerCell = 50000 
minCellsPergene = 5
for i, adata in enumerate(adatas):
    tid = tissueIdDict[str(i)]
    print("########################################")
    print("{0}) Processing {1}\n".format(i, tid))
    adata.var_names_make_unique()

    # Convert the sparse count matrices to dense represntation
    adata.X = adata.X.toarray()

    # Get individual plots dir for each sample
    indvBname    = tid
    indvPlotsDir = "{0}/01_preprocessing/{1}".format(plotsDir,tid); create_dir(indvPlotsDir)

    # 1.2.1) Calculate QC covariates
    print("- {0} with shape {1}".format(tid, adata.to_df().shape))
    sc.pp.calculate_qc_metrics(adata, inplace=True) # we now have many additional data types in the obs slot:
    adata.obs['n_counts']   = adata.X.sum(1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    adata.obs['n_genes']    = (adata.X > 0).sum(1)
    adata
    # obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'n_counts', 'log_counts', 'n_genes'
    # var: 'gene_ids', 'feature_types', 'genome', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'

    # Calculate mitochondrial/ribosomal genes fraction (percentage)
    # For each cell compute fraction of counts in mito/ribo genes vs. all genes
    mt_gene_mask         = [gene.startswith('mt-') for gene in adata.var_names]
    adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']
    rb_gene_mask         = [gene.startswith(("Rps","Rpl")) for gene in adata.var_names]
    adata.obs['rb_frac'] = adata.X[:, rb_gene_mask].sum(1)/adata.obs['n_counts']

    # 1.2.2) Plot QC metrics
    # Sample quality plots
    t1 = sc.pl.violin(adata, ['n_genes_by_counts', 'n_counts'], jitter=0.4, size=2, log=True, cut=0, show=False)
    plt.savefig("{0}/01_raw_{1}_nCounts_nGenes_plot.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')
    t2 = sc.pl.violin(adata, ['mt_frac','rb_frac'], jitter=0.4, size=2, log=False, cut=0, show=False)
    plt.savefig("{0}/01_raw_{1}_MtRbFrac_plot.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')

    # 1.2.3) Data quality summary plots
    p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', show=False)
    plt.savefig("{0}/01_raw_{1}_genes_counts_mtfrac_scatterplot.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')
    p2 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', show=False)
    plt.savefig("{0}/01_raw_{1}_genes_counts_mtfrac_scatterplot_zoomedin.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')

    p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='rb_frac', show=False)
    plt.savefig("{0}/01_raw_{1}_genes_counts_rbfrac_scatterplot.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')
    p2 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='rb_frac', show=False)
    plt.savefig("{0}/01_raw_{1}_genes_counts_rbfrac_scatterplot_zoomedin.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')

    # 1.2.4) Thresholdingecision based on counts
    p3 = sns.distplot(adata.obs['n_counts'], kde=False); #plt.show()
    plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')
    p4 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<2000], kde=False, bins=50); #plt.show()
    plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot_lessthan_2000.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')
    p5 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>5000], kde=False, bins=50); #plt.show()
    plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot_greaterthan_5000.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')

    # 1.2.5) Thresholding decision based on genes
    p6 = sns.distplot(adata.obs['n_genes'], kde=False, bins=50); # plt.show()
    plt.savefig("{0}/01_raw_{1}_genes_histogramplot.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')
    p7 = sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, bins=50); # plt.show()
    plt.savefig("{0}/01_raw_{1}_genes_histogramplot_lessthan_1000.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')

    # 1.2.6) Filter cells according to identified QC thresholds:
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

    # 1.2.7) Filter genes according to identified QC thresholds:
    print('Total number of genes: {:d}'.format(adata.n_vars))
    sc.pp.filter_genes(adata, min_cells=minCellsPergene)
    print('Number of genes after minCellsPergene filter: {:d}'.format(adata.n_vars))

    # 1.2.8) Compute variable genes
    # We first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.
    sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
    print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))

    # 1.2.9) Calculations for the visualizations
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
    sc.pp.neighbors(adata, random_state = 2105)
    sc.tl.umap(adata, random_state = 2105, n_components=3)

    # 1.2.10) Plot visualizations
    sc.pl.pca_scatter(adata, color='n_counts',show=False)
    plt.savefig("{0}/01_raw_{1}_clustering_ncounts_PCA.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')
    sc.pl.umap(adata, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
    plt.savefig("{0}/01_raw_{1}_clustering_UMAP.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')
    sc.pl.umap(adata, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
    plt.savefig("{0}/01_raw_{1}_clustering_UMAP_3D.png".format(indvPlotsDir, indvBname) , bbox_inches='tight', dpi=175); plt.close('all')

    # # 1.2.11) Expression recovery (denoising) the data on raw counts
    # sce.pp.dca(adata)

    # 1.2.12) Add to the filtered adatas list
    fadatas.append(adata)

# 1.3) Merge 10x datasets for different mices
# https://github.com/theislab/scanpy/issues/267
adata = fadatas[0].concatenate(fadatas[1:])

# 1.4) Make variable names unique
adata.var_names_make_unique()

# 1.5) Add tissue id column for the batches
adata.obs['tissueID'] = adata.obs['batch'].map(tissueIdDict)

condition_key = "tissueID"
cell_type_key = "cellType"

# 1.6) Checking the total size of the data set
adata.shape # We have 5887 cells and 11029 genes in the dataset

# Reset the parameters for merged data
minGenesPerCell = 150
minCountPerCell = 100
maxCountPerCell = 50000 
minCellsPergene = 75

# 2) Plot QC matrices for raw filtered merged data
# 2.1) Plot sample QC metrics for merged data
t1 = sc.pl.violin(adata, 'n_counts', groupby='tissueID', size=2, log=True, cut=0, show=False)
plt.savefig("{0}/01_raw_{1}_tissueID_nCounts_plot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
t2 = sc.pl.violin(adata, 'mt_frac', groupby='tissueID', show=False)
plt.savefig("{0}/01_raw_{1}_tissueID_mtFraction_plot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.2) Data quality summary plots
p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/01_raw_{1}_genes_counts_scatterplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p2 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', show=False)
plt.savefig("{0}/01_raw_{1}_genes_counts_scatterplot_zoomedin.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.3) QC plots based on counts
p3 = sns.distplot(adata.obs['n_counts'], kde=False); #plt.show()
plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p4 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<2000], kde=False, bins=50); #plt.show()
plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot_lessthan_2000.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p5 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>5000], kde=False, bins=50); #plt.show()
plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot_greaterthan_5000.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.4) QC plots based on genes
p6 = sns.distplot(adata.obs['n_genes'], kde=False, bins=50); # plt.show()
plt.savefig("{0}/01_raw_{1}_genes_histogramplot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
p7 = sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, bins=50); # plt.show()
plt.savefig("{0}/01_raw_{1}_genes_histogramplot_lessthan_1000.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 2.5) Filter cells according to identified QC thresholds:
origadata = adata.copy()
print('Total number of cells: {:d}'.format(adata.n_obs))
sc.pp.filter_cells(adata, min_counts = minCountPerCell)
print('Number of cells after min count filter: {:d}'.format(adata.n_obs))
sc.pp.filter_cells(adata, max_counts = 40000)
print('Number of cells after max count filter: {:d}'.format(adata.n_obs))
sc.pp.filter_cells(adata, min_genes = minGenesPerCell)
print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

# Total number of cells: 5887
# Number of cells after min count filter: 5887
# Number of cells after max count filter: 5885
# Number of cells after gene filter: 5842

# 2.6) Filter genes according to identified QC thresholds:
# Min minCellsPergene cells - filters out 0 count genes
print('Total number of genes: {:d}'.format(adata.n_vars))
sc.pp.filter_genes(adata, min_cells=minCellsPergene)
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
# Total number of genes: 11029
# Number of genes after cell filter: 9606

# 2.7) Compute highly variable genes (HVGs)
# Next, we first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000, batch_key = 'tissueID')
print("Highly variable genes intersection: %d"%sum(adata.var.highly_variable_intersection))
print("Number of batches where gene is variable:")
print(adata.var.highly_variable_nbatches.value_counts())
# Highly variable genes intersection: 1096
# Number of batches where gene is variable:
# 1    2921
# 2    2209
# 0    1955
# 3    1425
# 4    1096
# Name: highly_variable_nbatches, dtype: int64

# Calculations for the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(adata, random_state = 2105)
sc.tl.umap(adata, random_state = 2105, n_components=3)

# 2.9) Plot visualizations
# PCA
sc.pl.pca_scatter(adata, color='tissueID', components = ['1,2','2,3','3,4','4,5','5,6','6,7'], ncols=3, hspace=0.35, wspace=0.35, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/01_raw_{1}_tissueID_PCA.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# UMAPS
fig = plt.figure(figsize=(16,6))
fig.suptitle('Raw data')
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  
sc.pl.umap(adata, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(1, 2, 2, projection='3d'); 
sc.pl.umap(adata, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/01_raw_{1}_tissueID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 5.8) Save the trVAE batch corrected adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata  = sc.read_h5ad(adatafile)

########################
rawadata = adata.copy()
# 3) # Normalization using SCRAN
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
size_factors = computeSumFactors(data_mat, clusters=input_groups, min.mean=0.25)
# Delete adata_pp
del adata_pp
# Visualize the estimated size factors
adata.obs['size_factors'] = size_factors
fig = plt.figure(figsize=(16,6))
fig.suptitle('Estimated size factors')
ax = fig.add_subplot(1, 2, 1)
sc.pl.scatter(adata, 'size_factors', 'n_counts', ax=ax, show=False)
ax = fig.add_subplot(1, 2, 2)
sc.pl.scatter(adata, 'size_factors', 'n_genes', ax=ax, show=False)
plt.tight_layout()
plt.savefig("{0}/02_norm_{1}_scran_sizefactors_plots.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
# Keep the count data in a counts layer
adata.layers["counts"] = adata.X.copy()
# Normalize adata 
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)
# Store the full data set in 'raw' as log-normalised data for statistical testing
adata.raw = adata
#########################################################################
# 4) Biological correction
# 4.1) Read cell cycle genes
cc_genes         = pd.read_table(ccGenes_macosko, delimiter='\t')
s_genes          = cc_genes['S'].dropna()
g2m_genes        = cc_genes['G2.M'].dropna()
# For mouse only
s_genes_mm       = [gene.lower().capitalize() for gene in s_genes]
g2m_genes_mm     = [gene.lower().capitalize() for gene in g2m_genes]
s_genes_mm_ens   = adata.var_names[np.in1d(adata.var_names, s_genes_mm)]
g2m_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, g2m_genes_mm)]
# 4.2) Score cell cycle genes
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)

# 4.3) Visualize the effects of cell cycle
fig = plt.figure(figsize=(16,12))
fig.suptitle('Effects of Cell Cycle')
ax = fig.add_subplot(2, 2, 1)
sc.pl.umap(adata, color=['S_score']  , ax=ax, use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
ax = fig.add_subplot(2, 2, 2)
sc.pl.umap(adata, color=['G2M_score'], ax=ax, use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
ax = fig.add_subplot(2, 2, 3)
sc.pl.umap(adata, color='phase', ax=ax, use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
ax = fig.add_subplot(2, 2, 4, projection='3d')
sc.pl.umap(adata, color='phase', ax=ax, use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.tight_layout()
plt.savefig("{0}/02_norm_{1}_scran_cell_cycle_plots.png".format(plotsDir, bname) , bbox_inches='tight', dpi=750); plt.close('all')

# 4.4) Save the normalized cell cycle corrected adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/02_normCC_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/02_normCC_{1}_adata.h5ad" .format(dataDir, projName); normadata  = sc.read_h5ad(adatafile)
#########################################################################
normadata = adata.copy()
rawadatas = adatas.copy()

# adata  = normadata.copy()
# adatas = rawadata.copy()
# 5) Technical correction: Batch Correction using TrVAE

# 5.1) Normalizing & Extracting Top 1000 Highly Variable Genes
# We can preserve more genes (up to 7000 like scGen) but in order to train the network quickly, we will extract top 1000 genes
adata2 = adata.copy()
sc.pp.log1p(adata2)
sc.pp.highly_variable_genes(adata2, n_top_genes=1000)
adata2 = adata2[:, adata2.var['highly_variable']]
adata2

# 5.2) Train/Test Split
train_adata, valid_adata = trvae.utils.train_test_split(adata2, train_frac=0.80)
train_adata.shape, valid_adata.shape
# ((4348, 7109), (1087, 7109))

# 5.3) Calculate number of batches
n_conditions = len(train_adata.obs[condition_key].unique().tolist())
conditions = adata2.obs[condition_key].unique().tolist()
condition_encoder = trvae.utils.create_dictionary(conditions, [])
condition_encoder
# {'bulk1001': 0, 'bulk997': 1, 'bulk1018': 2, 'stomach1001': 3}

# 5.4)Create the network
# Some of network parameters:
# # x_dimension: number of features (necessary)
# # n_conditons: total number of batches (necessary)
# # architecture: architecture of the network (optional)
# # output_activation: activation function of trVAE's last layer
# # alpha: coefficient of KL Divergence loss (optional)
# # beta: coefficient of MMD loss (optional)
# # eta: coefficient of reconstruction (MSE or SSE) loss (optional) can be one of the relu, leaky_relu, linear, ...
network = trvae.archs.trVAE(x_dimension=train_adata.shape[1],
                            architecture=[128, 32],
                            z_dimension=10,
                            n_conditions=n_conditions,
                            alpha=0.00005,
                            beta=50,
                            eta=100,
                            output_activation='relu')

# 5.5) Training trVAE
# Parameters:
# train_adata: anndata object for training trVAE
# valid_adata: anndata object for validating trVAE (default = None)
# condition_key: (str) name of the column containing batches' names in train_adata and valid_adata (Necessary)
# condition_encoder: (dict) dictionary of encoded batches (keys: batch names, values: encoded integers) (default = None)
# verbose: (int) level of verbosity (default = 0)
# NOTE: This will take 8 to 10 minutes on WS4 with single CPU
outModelFile    = "{0}/trVAE_{1}_network.model".format(dataDir,projName)
network.train(train_adata,
              valid_adata,
              condition_encoder,
              condition_key,
              n_epochs=1000,
              batch_size=1024,
              verbose=5,
              early_stop_limit=750,
              lr_reducer=0,
              shuffle=True,
              save=True,
              )

# 5.6) Getting batch-corrected adata
# Now two matrices have been added to adata
# mmd_latent: (numpy ndarray) output of MMD Layer in trVAE
# reconstructed: (numpy ndarray) reconstructed data with dimension of original feature space
# z_latent: (numpy ndarray) output of bottleneck layer of trVAE (optional)
# For evaluating how good trVAE has corrected the batches, we recommend using mmd_latent matrix.
labels, _ = trvae.tl.label_encoder(adata2, condition_key=condition_key, label_encoder=condition_encoder)
network.get_corrected(adata2, labels, return_z=True)

# 5.7) MMD Layer UMAP visualization
# mmd_latent = adata2.obsm['mmd_latent']
# mmd_adata = sc.AnnData(mmd_latent, obs=adata2.obs)
# mmd_adata
adata.obsm['mmd_latent']    = adata2.obsm['mmd_latent']
adata.obsm['z_latent']      = adata2.obsm['z_latent']
adata.obsm['reconstructed'] = adata2.obsm['reconstructed']
# sc.pp.neighbors(adata, random_state = 2105, n_neighbors=10, use_rep = "mmd_latent")
sc.pp.neighbors(adata, random_state = 2105, n_neighbors=7, use_rep = "mmd_latent")
sc.tl.umap(adata, random_state = 2105, n_components=3)
fig = plt.figure(figsize=(16,13))
fig.suptitle('TissueID')
# 2D projection
ax = fig.add_subplot(2, 2, 1);                  sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 2, 2);                  sc.pl.umap(adata, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE UMAP")
# 3D projection
ax = fig.add_subplot(2, 2, 3, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 2, 4, projection='3d'); sc.pl.umap(adata                , ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE UMAP")
plt.tight_layout()
plt.savefig("{0}/03_norm_TrVAE_batchCorrection_{1}_tissueID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 5.8) Save the trVAE batch corrected adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/03_norm_TrVAE_batchCorrection_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/03_norm_TrVAE_batchCorrection_{1}_adata.h5ad" .format(dataDir, projName); trvaeBCadata  = sc.read_h5ad(adatafile)

############################################################
############################################################
# Using TrVAE for batch correction
############################################################
trvaeBCadata = adata.copy() # (5885, 11023)
# adata = trvaeBCadata.copy() # (5885, 11023)

# 6) Plot UMAP for all the tumor cells 
# Color the cells that have human myc and ires
cellBarCodes = pd.read_csv('/media/rad/HDD2/temp_manec/hgMycIresCd2_cellIDs.txt', sep="\t", header=None).values.tolist()
cl  = sum(cellBarCodes, [])
ucl = get_unique_list(cl)
# In [34]: len(ucl)
# Out[34]: 1832

mylist = adata.obs.index.values
humaniresmyc = list()
for e in mylist: 
  flag = 0
  for s in ucl: 
      if s in e: 
          flag = 1 
          break
  humaniresmyc.append(flag)

adata.obs['hgMycIresCd2'] = humaniresmyc
fig = plt.figure(figsize=(16,6))
fig.suptitle('hgMycIresCd2')
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  
sc.pl.umap(adata, legend_loc=None, ax=ax, color="hgMycIresCd2", color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(1, 2, 2, projection='3d'); 
sc.pl.umap(adata, ax=ax, color="hgMycIresCd2", color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/03_normTrVAE_{1}_Tumor_hgMycIresCd2_CellIDs_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# --------------------------------------------------------
# Get tumors for all individual vectors ('humanMyc', 'gap', 'ires', 'humanCd2')
# Unique cell barcode list of list
ucblofl = list()
ucbd    = defaultdict(list)
for t in ['humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2']:
  # Cell barcode data frame
  cbDF = pd.read_csv('/media/rad/HDD2/temp_manec/hgMycIresCd2_{0}_cellIDs.txt'.format(t), sep="\t", header=None).values.tolist()
  # Unique cell barcode list
  ucbl = get_unique_list(sum(cbDF, []))
  ucbd[t] = ucbl
  ucblofl.append(ucbl)
  print(len(ucbl))

  # Add the flag to adata (Unique cell barcode flag list)
  ucbfl = list()
  for e in mylist: 
    flag = 0
    for s in ucbl: 
        if s in e: 
            flag = 1 
            break
    ucbfl.append(flag)
  adata.obs[t] = ucbfl

import upsetplot
from upsetplot import from_memberships, from_contents
upsetContent = from_contents(ucbd)
memberships  = from_memberships(ucblofl)

fig = plt.figure(figsize=(30,10))
upsetplot.plot(upsetContent, sum_over=False,show_counts='%d', fig=fig)
plt.savefig("{0}/03_normTrVAE_{1}_Tumor_hgMycIresCd2_individual_CellIDs_upsetPlot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Plot all hgMycIresCd2 vector components in separate UMAP
fig = plt.figure(figsize=(16,6))
fig.suptitle('hgMycIresCd2 vector components')
# 2D projection
sc.pl.umap(adata, color=['hgMycIresCd2','humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2'], color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, show=False)
plt.savefig("{0}/03_normTrVAE_{1}_Tumor_hgMycIresCd2_Components_UMAP_2D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
# 3D projection
sc.pl.umap(adata, color=['hgMycIresCd2','humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2'], color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/03_normTrVAE_{1}_Tumor_hgMycIresCd2_Components_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Get tumor cells per tissue (and in percentage)
tumorDF = adata.obs[['batch','humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2','hgMycIresCd2']].copy()
# Replace batch number with batch names
tumorDF.replace({'batch': tissueIdDict}, inplace=True)
# Remove index for groupby
tumorDF.reset_index(drop=True, inplace=True)

# Get the number of cells for each cluster in every tissue
ncellstumorDF = tumorDF.groupby(['batch']).sum()
ncellstumorDF.to_csv("{0}/03_{1}_tumorComponents_per_tissueid.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="tissueID")




##################################################################
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

# Number of cells in each cluster
# adata.obs['louvain_r1.5'].value_counts()                                                                                                                     # 0     821

# Plot visualizations
# Visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_normTrVAE_{1}_clustering_all_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/03_normTrVAE_{1}_clustering_all_louvain_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

cluster_key   = "louvain_r0.5"
cluster_bname = "louvain_r05"
fig = plt.figure(figsize=(32,8))
# 2D projection
ax = fig.add_subplot(2, 5, 1);                  sc.pl.umap(rawadata,                  ax=ax, color="tissueID"  , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw tissueID UMAP")
ax = fig.add_subplot(2, 5, 2);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="tissueID"  , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE tissueID UMAP")
ax = fig.add_subplot(2, 5, 3);                  sc.pl.umap(adata   ,                  ax=ax, color=cluster_key   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format(cluster_key))
ax = fig.add_subplot(2, 5, 4);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="log_counts", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="log_counts UMAP")
ax = fig.add_subplot(2, 5, 5);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="mt_frac UMAP")
# 3D projection
ax = fig.add_subplot(2, 5, 6, projection='3d'); sc.pl.umap(rawadata, legend_loc=None,  ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw tissueID UMAP")
ax = fig.add_subplot(2, 5, 7, projection='3d'); sc.pl.umap(adata   , legend_loc=None,  ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE tissueID UMAP")
ax = fig.add_subplot(2, 5, 8, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format(cluster_key))
ax = fig.add_subplot(2, 5, 9, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color="log_counts"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="log_counts UMAP")
ax = fig.add_subplot(2, 5, 10, projection='3d'); sc.pl.umap(adata  , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="mt_frac UMAP")
plt.tight_layout()
plt.savefig("{0}/03_normTrVAE_{1}_{2}_tissueID_counts_mtfrac_UMAP.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Louvain UMAPs
fig = plt.figure(figsize=(16,6))
fig.suptitle("{0} UMAP".format(cluster_key))
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  
sc.pl.umap(adata, legend_loc=None, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(1, 2, 2, projection='3d'); 
sc.pl.umap(adata, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/03_normTrVAE_{1}_clustering_{2}_UMAP_2D3D.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
cluster_key_groups = adata.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = adata.obs[cluster_key].value_counts().to_dict()

# Louvain UMAPs
ncols  = len(cluster_key_groups) + 1
fig = plt.figure(figsize=(7*ncols, 14))
fig.suptitle("{0} UMAP".format(cluster_key))
# Main Louvain Cluster
ax = fig.add_subplot(2, ncols, 1); sc.pl.umap(adata, legend_loc=None, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(2, ncols, 2, projection='3d'); sc.pl.umap(adata, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
# Partial visualizaton of a subset of groups in embedding
m=3; n=4
for i,b in enumerate(cluster_key_groups):
  print(i, b)
  ax = fig.add_subplot(2, ncols, i+m);                  sc.pl.umap(adata, legend_loc=None, ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
  ax = fig.add_subplot(2, ncols, i+n, projection='3d'); sc.pl.umap(adata                 , ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
  m+=1; n+=1

plt.tight_layout()
plt.savefig("{0}/03_normTrVAE_{1}_clustering_{2}_UMAP_individual_clusters.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 7.2) Plot separate bar plots, coloured in by cluster annotation, for each tissue
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
plt.savefig("{0}/03_normTrVAE_{1}_clustering_{2}_tissueID_cluster_barplot.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')


##################################################################
# 8.1) Marker genes & cluster annotation
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby=cluster_key, key_added='rank_genes')

# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes', fontsize=12, show=False)
plt.savefig("{0}/03_normTrVAE_{1}_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# Annotation of cluster r_0.5 with known marker genes
markerDir = "{0}/markerDir".format(plotsDir); create_dir(markerDir)
subplot_title_fontsize = 12
subplot_title_width    = 50

# Read the marker genes into a pandas dataframe
marker_file  = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_V2.txt'
markersDF    = pd.read_csv(marker_file, sep="\t")
marker_genes = markersDF.groupby('CellLines')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()
marker_genes_cellTypes = markersDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()

# markerexp = marker_gene_expression(adata, marker_genes_cellTypes, gene_symbol_key=None, partition_key='louvain_r0.5')

# For mouse cell atlas marker genes
ma_marker_file       = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_mouse_cellatlas_V1.txt'
ma_markersDF         = pd.read_csv(ma_marker_file, sep="\t", header=None, index_col=None)
ma_markersDF         = ma_markersDF[0].str.split(",", n = 1, expand = True)
ma_markersDF.columns = ['CellTypes', 'MarkerGenes']
ma_marker_genes      = ma_markersDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()

# l = ['a', 'b', 'c', 'd', 'f']
# d = {'A':['a','x','c'], 'B':['c','d'],'C':['x','y']}

# for k,v in d.items():
#   nv = [x for x in v if x in l]
#   d[k] = nv

# Get all the gene names in the adata object
genespresent = adata.var.index.values.tolist()

# Generate the UMAPs for each marker categories
for k,v in marker_genes_cellTypes.items():
  print("\n- Original list {0}: {1}".format(k,v))
  validgenes = [x for x in v if x in genespresent]
  ids = np.in1d(adata.var_names,validgenes)
  print("- Genes present {0}: {1}".format(k,validgenes))

  ngenes = len(validgenes)
  nrows  = ngenes + 2
  adata.obs['{0}_marker_expr'.format(k)] = adata.X[:,ids].mean(1)

  fig = plt.figure(figsize=(14,6*nrows))
  # fig.suptitle('Stomach_marker_list_V1')
  # Plot cluster
  ax = fig.add_subplot(nrows, 2, 1);                  sc.pl.umap(adata, legend_loc='on data', ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format(cluster_key))
  ax = fig.add_subplot(nrows, 2, 2, projection='3d'); sc.pl.umap(adata                      , ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format(cluster_key))

  # Plots mean marker genes
  ax = fig.add_subplot(nrows, 2, 3);                  sc.pl.umap(adata, legend_loc=None     , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))
  ax = fig.add_subplot(nrows, 2, 4, projection='3d'); sc.pl.umap(adata                      , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))

  # Plot individual marker genes
  m=5; n=6
  for i,mgene in enumerate(validgenes):
    # print(i+m, i+n, mgene)
    ax = fig.add_subplot(nrows, 2, i+m);                  sc.pl.umap(adata, legend_loc=None     , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows, 2, i+n, projection='3d'); sc.pl.umap(adata                      , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
    m+=1; n+=1

  plt.tight_layout()
  plt.savefig("{0}/21_{1}_marker_genes_stomach_V2_{2}_UMAPs.png".format(markerDir, bname, k) , bbox_inches='tight', dpi=100); plt.close('all')

# Generate the UMAPs for each marker categories
for k,v in ma_marker_genes.items():
  print("\n- Original list {0}: {1}".format(k,v))
  validgenes = [x for x in v if x in genespresent]
  ids = np.in1d(adata.var_names,validgenes)
  print("- Genes present {0}: {1}".format(k,validgenes))

  ngenes = len(validgenes)
  nrows  = ngenes + 2
  adata.obs['{0}_marker_expr'.format(k)] = adata.X[:,ids].mean(1)

  fig = plt.figure(figsize=(14,6*nrows))
  # fig.suptitle('Stomach_marker_list_V1')
  # Plot cluster
  ax = fig.add_subplot(nrows, 2, 1);                  sc.pl.umap(adata, legend_loc='on data', ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format(cluster_key))
  ax = fig.add_subplot(nrows, 2, 2, projection='3d'); sc.pl.umap(adata                      , ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format(cluster_key))

  # Plots mean marker genes
  ax = fig.add_subplot(nrows, 2, 3);                  sc.pl.umap(adata, legend_loc=None     , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))
  ax = fig.add_subplot(nrows, 2, 4, projection='3d'); sc.pl.umap(adata                      , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))

  # Plot individual marker genes
  m=5; n=6
  for i,mgene in enumerate(validgenes):
    # print(i+m, i+n, mgene)
    ax = fig.add_subplot(nrows, 2, i+m);                  sc.pl.umap(adata, legend_loc=None     , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows, 2, i+n, projection='3d'); sc.pl.umap(adata                      , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
    m+=1; n+=1

  plt.tight_layout()
  plt.savefig("{0}/32_{1}_mouse_cellatlas_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, k) , bbox_inches='tight', dpi=100); plt.close('all')

# Save the louvain information in external file
louvainsDF = pd.DataFrame(adata.obs[cluster_key])
louvainsDF.to_csv("{0}/03_{1}_louvains.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="cellId")

# 7.5) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); markeradata  = sc.read_h5ad(adatafile)

################################################################################
# adata = markeradata.copy()

# 8) CellType assignments
# 8.1) Add new categories 
# Categories to rename
adata.obs[cluster_key].cat.categories
# Get a new cell type column from the annotation of the louvain_r0.5 clusters
adata.obs['cellType'] = adata.obs[cluster_key]
# Add new categories
adata.obs['cellType'].cat.add_categories(['Dendritic_Macrophages','Unknown1','Unknown2','Erythrocytes','Endothelial_Epithelial_Igfbp3pos','Fibroblasts','ECL_Restin_Macrophages','Progenitor_at_neck','Tumor'], inplace=True) 
# Get a new subcluster column
# 0           = 'Dendritic_Macrophages'
# 1           = 'Unknown0'
# 3           = 'Unknown1'
# 4           = 'Erythrocytes'
# 5           = 'Endothelial/Epithelial_Igfbp3⁺'
# 7           = 'Fibroblasts'
# 8           = 'ECL_Restin_Macrophages'
# 10          = 'Progenitor_at_neck'
# 2,6,9,11    = 'Tumor'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='0' ]  = 'Dendritic_Macrophages'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='1' ]  = 'Unknown1'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='3' ]  = 'Unknown2'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='4' ]  = 'Erythrocytes'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='5' ]  = 'Endothelial_Epithelial_Igfbp3pos'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='7' ]  = 'Fibroblasts'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='8' ]  = 'ECL_Restin_Macrophages'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='10']  = 'Progenitor_at_neck'
adata.obs['cellType'].loc[(adata.obs['cellType']=='2')|(adata.obs['cellType']=='6')|(adata.obs['cellType']=='9')|(adata.obs['cellType']=='11')]  = 'Tumor'
# Remove old categories
adata.obs['cellType'].cat.remove_categories(adata.obs[cluster_key].cat.categories.tolist(), inplace=True)
# List new categories
adata.obs['cellType'].cat.categories
# Draw Umaps with the categories
fig = plt.figure(figsize=(36,20))
fig.suptitle('CellType UMAP')
# 2D projection
ax = fig.add_subplot(2, 3, 1);                  sc.pl.umap(adata, legend_loc='on data', ax=ax, color="cellType", palette=sc.pl.palettes.vega_20, size=150, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(2, 3, 2);                  sc.pl.umap(adata, legend_loc=None     , ax=ax, color="cellType".format(k),  palette=sc.pl.palettes.vega_20, size=150, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(2, 3, 3, projection='3d'); sc.pl.umap(adata                      , ax=ax, color="cellType", palette=sc.pl.palettes.vega_20, size=150, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
# Save the UMAP
plt.savefig("{0}/04_normTrVAE_{1}_clustering_CellType_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=200); plt.close('all')

# 8.2) Plot separate bar plots, coloured in by cluster annotation, for each tissue
# Convert palette into colormap
clcmap = ListedColormap(sc.pl.palettes.vega_20)
# Get the DF of tissue and clusters
clusterBatchDF = adata.obs[['batch','cellType']].copy()
# Replace batch number with batch names
clusterBatchDF.replace({'batch': tissueIdDict}, inplace=True)
# Remove index for groupby
clusterBatchDF.reset_index(drop=True, inplace=True)
# Get the number of cells for each cluster in every tissue
ncellsClusterBatchDF = clusterBatchDF.groupby(['batch','cellType']).size()
# Get the percent of cells for each cluster in every tissue 
pcellsClusterBatchDF = pd.crosstab(index=clusterBatchDF['batch'], columns=clusterBatchDF['cellType'], values=clusterBatchDF['cellType'], aggfunc='count', normalize='index')

# Plot the barplots
fig = plt.figure(figsize=(16,6)); fig.suptitle("Cells for each {0} in each tissue".format('cellType'))
# plot numbers of cells
ax = fig.add_subplot(1, 2, 1); ncellsClusterBatchDF.unstack().plot(kind='barh', stacked=True, colormap=clcmap, ax=ax, legend=None, title="Number of cells")
# plot percent of cells
ax = fig.add_subplot(1, 2, 2); pcellsClusterBatchDF.plot(kind='barh',stacked=True, colormap=clcmap, ax=ax, title="% of cells")
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='cellType', title_fontsize=12)
plt.tight_layout() # For non-overlaping subplots
plt.savefig("{0}/04_subGroup_{1}_clustering_{2}_tissueID_cellType_barplot.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 8.3) Calculate marker genes (one vs. rest)
sc.tl.rank_genes_groups(adata, groupby='cellType', key_added='rank_genes')
# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes', fontsize=12, show=False)
plt.savefig("{0}/04_subGroup_{1}_{2}_marker_genes_ranking_cellType.png".format(plotsDir, bname, 'cellType') , bbox_inches='tight', dpi=175); plt.close('all')

# 8.4) Calculate pairwise marker genes list
# Get all cellTypes into the list
cellTypeCategories = adata.obs['cellType'].cat.categories.tolist()
# Get the list of all unique pairwise combinations
cellTypePairComb = [comb for comb in combinations(cellTypeCategories, 2)]
# Get pairwise plots dir
pwPlotsDir = "{0}/rankedGenes/allCategories".format(plotsDir); create_dir(pwPlotsDir)
# Get pairwise data dir
pwDataDir  = "{0}/rankedGenes/allCategories".format(dataDir); create_dir(pwDataDir)
# Calculate pairwise marker genes  
for grp,ref in cellTypePairComb:
  print("- Calculating pairwise marker genes for group_v_reference: {0}_v_{1}".format(grp, ref))
  # Get genes ranking
  keyName = 'rank_genes_{0}_v_{1}'.format(grp,ref)
  sc.tl.rank_genes_groups(adata, groupby='cellType', groups=[grp], key_added=keyName, reference=ref, n_genes=adata.shape[1])
  # Plot top 20 ranked genes
  sc.pl.rank_genes_groups(adata, key=keyName, groups=[grp], fontsize=12, show=False)
  # Save it in a figure
  plt.savefig("{0}/04_{1}_all_cellType_{2}_v_{3}.png".format(pwPlotsDir, bname, grp, ref) , bbox_inches='tight'); plt.close('all')
  # Get the dataframe of DE parameters
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(adata.uns[keyName][n])[grp]
  # Add treatment and reference group name
  ngDF['Treatment'] = grp
  ngDF['Reference'] = ref
  # Save the dataframe
  ngDF.to_csv("{0}/04_{1}_{2}.txt".format(pwDataDir, projName, keyName), sep='\t', header=True, index=False, float_format='%.2g')

# 8.5) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/05_cellType_assigned_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)

# ################################# 
#     SubCluster Analysis 
# #################################

# # Read back the corrected adata object
# adatafile  = "{0}/05_cellType_assigned_{1}_adata.h5ad" .format(dataDir, projName); cellTypeadata  = sc.read_h5ad(adatafile)
# adata = cellTypeadata.copy()

# 9.1) Get the sub group ann data object 
# Remove the cells that are not needed
# 2           = 'Dendritic'
# 3           = 'Endothelial_Epithelial_Igfbp3pos'
# 5           = 'Erythrocytes'
# 7           = 'Macrophages'
# 9           = 'Fibroblasts'
adataSubGroup = adata[~((adata.obs['cellType']=='Dendritic_Macrophages') | (adata.obs['cellType']=='Endothelial_Epithelial_Igfbp3pos') |(adata.obs['cellType']=='Erythrocytes') | (adata.obs['cellType']=='Fibroblasts'))].copy()
# adataSubGroup.shape # (3483, 9606)

# Calculations for the visualizations
sc.pp.neighbors(adataSubGroup, random_state = 2105, n_neighbors=7, use_rep = "mmd_latent")
sc.tl.umap(adataSubGroup, random_state = 2105, n_components=3)
fig = plt.figure(figsize=(16,6))
fig.suptitle('TissueID')
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE UMAP")
# 3D projection
ax = fig.add_subplot(1, 2, 2, projection='3d'); sc.pl.umap(adataSubGroup                , ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE UMAP")
plt.tight_layout()
plt.savefig("{0}/05_subGroup_{1}_tissueID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 9.2) Perform clustering - using highly variable genes
sc.tl.louvain(adataSubGroup, key_added='louvain', random_state=2105)
sc.tl.louvain(adataSubGroup, resolution=1, key_added='louvain_r1', random_state=2105)
sc.tl.louvain(adataSubGroup, resolution=1.5, key_added='louvain_r1.5', random_state=2105)
sc.tl.louvain(adataSubGroup, resolution=2.0, key_added='louvain_r2', random_state=2105)
for i in np.linspace(0.1,0.9,9):
    try:
        sc.tl.louvain(adataSubGroup, resolution=i, key_added='louvain_r{0}'.format(i), random_state=2105)
        print(adataSubGroup.obs['louvain_r{0:0.1f}'.format(i)].value_counts())
    except:
        print("- Error in r: {0}".format(i))
sc.tl.louvain(adataSubGroup, resolution=0.3, key_added='louvain_r0.3', random_state=2105)
sc.tl.louvain(adataSubGroup, resolution=0.7, key_added='louvain_r0.7', random_state=2105)

# Visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adataSubGroup, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/05_subGroup_{1}_clustering_all_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adataSubGroup, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/05_subGroup_{1}_clustering_all_louvain_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Plot tumor cells
fig = plt.figure(figsize=(16,6))
fig.suptitle('hgMycIresCd2')
ax = fig.add_subplot(1, 2, 1); sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color="hgMycIresCd2", color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(1, 2, 2, projection='3d'); sc.pl.umap(adataSubGroup, ax=ax, color="hgMycIresCd2", color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/05_subGroup_{1}_Tumor_hgMycIresCd2_CellIDs_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Subcluster keys
cluster_key            = "louvain_r0.5"
cluster_bname          = "louvain_r05"
subplot_title_fontsize = 12
subplot_title_width    = 50

# Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
cluster_key_groups = adataSubGroup.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = adataSubGroup.obs[cluster_key].value_counts().to_dict()

# Louvain UMAPs
ncols  = len(cluster_key_groups) + 1
fig = plt.figure(figsize=(7*ncols, 14))
fig.suptitle("{0} UMAP".format(cluster_key))
# Main Louvain Cluster
ax = fig.add_subplot(2, ncols, 1); sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(2, ncols, 2, projection='3d'); sc.pl.umap(adataSubGroup, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
# Partial visualizaton of a subset of groups in embedding
m=3; n=4
for i,b in enumerate(cluster_key_groups):
  print(i, b)
  ax = fig.add_subplot(2, ncols, i+m);                  sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
  ax = fig.add_subplot(2, ncols, i+n, projection='3d'); sc.pl.umap(adataSubGroup                 , ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
  m+=1; n+=1

plt.tight_layout()
plt.savefig("{0}/05_subGroup_{1}_clustering_{2}_UMAP.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 9.3) Plot separate bar plots, coloured in by cluster annotation, for each tissue
# Convert palette into colormap
clcmap = ListedColormap(sc.pl.palettes.vega_20)
# Get the DF of tissue and clusters
clusterBatchDF = adataSubGroup.obs[['batch','{0}'.format(cluster_key)]].copy()
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
plt.savefig("{0}/05_subGroup_{1}_clustering_{2}_tissueID_cluster_barplot.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 9.4) Marker genes & cluster annotation
# Calculate marker genes
sc.tl.rank_genes_groups(adataSubGroup, groupby=cluster_key, key_added='rank_genes')

# Plot marker genes
sc.pl.rank_genes_groups(adataSubGroup, key='rank_genes', fontsize=12, show=False)
plt.savefig("{0}/05_subGroup_{1}_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# Annotation of cluster r_0.5 with known marker genes
markerDir = "{0}/subGroup/markerDir".format(plotsDir); create_dir(markerDir)
subplot_title_fontsize = 12
subplot_title_width    = 50

# Read the marker genes into a pandas dataframe
marker_file  = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_V2.txt'
markersDF    = pd.read_csv(marker_file, sep="\t")
marker_genes = markersDF.groupby('CellLines')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()
marker_genes_cellTypes = markersDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()

# For mouse cell atlas marker genes
ma_marker_file       = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_mouse_cellatlas_V1.txt'
ma_markersDF         = pd.read_csv(ma_marker_file, sep="\t", header=None, index_col=None)
ma_markersDF         = ma_markersDF[0].str.split(",", n = 1, expand = True)
ma_markersDF.columns = ['CellTypes', 'MarkerGenes']
ma_marker_genes      = ma_markersDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()

# Get all the gene names in the adataSubGroup object
genespresent = adataSubGroup.var.index.values.tolist()

# Generate the UMAPs for each marker categorie
for k,v in marker_genes_cellTypes.items():
  print("\n- Original list {0}: {1}".format(k,v))
  validgenes = [x for x in v if x in genespresent]
  ids = np.in1d(adataSubGroup.var_names,validgenes)
  print("- Genes present {0}: {1}".format(k,validgenes))

  ngenes = len(validgenes)
  nrows  = ngenes + 2
  adataSubGroup.obs['{0}_marker_expr'.format(k)] = adataSubGroup.X[:,ids].mean(1)

  fig = plt.figure(figsize=(14,6*nrows))
  # fig.suptitle('Stomach_marker_list_V1')
  # Plot cluster
  ax = fig.add_subplot(nrows, 2, 1);                  sc.pl.umap(adataSubGroup, legend_loc='on data', ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format(cluster_key))
  ax = fig.add_subplot(nrows, 2, 2, projection='3d'); sc.pl.umap(adataSubGroup                      , ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format(cluster_key))

  # Plots mean marker genes
  ax = fig.add_subplot(nrows, 2, 3);                  sc.pl.umap(adataSubGroup, legend_loc=None     , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))
  ax = fig.add_subplot(nrows, 2, 4, projection='3d'); sc.pl.umap(adataSubGroup                      , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))

  # Plot individual marker genes
  m=5; n=6
  for i,mgene in enumerate(validgenes):
    # print(i+m, i+n, mgene)
    ax = fig.add_subplot(nrows, 2, i+m);                  sc.pl.umap(adataSubGroup, legend_loc=None     , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows, 2, i+n, projection='3d'); sc.pl.umap(adataSubGroup                      , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
    m+=1; n+=1

  plt.tight_layout()
  plt.savefig("{0}/21_{1}_marker_genes_stomach_V2_{2}_UMAPs.png".format(markerDir, bname, k) , bbox_inches='tight', dpi=100); plt.close('all')

# Generate the UMAPs for each marker categories
for k,v in ma_marker_genes.items():
  print("\n- Original list {0}: {1}".format(k,v))
  validgenes = [x for x in v if x in genespresent]
  ids = np.in1d(adataSubGroup.var_names,validgenes)
  print("- Genes present {0}: {1}".format(k,validgenes))

  ngenes = len(validgenes)
  nrows  = ngenes + 2
  adataSubGroup.obs['{0}_marker_expr'.format(k)] = adataSubGroup.X[:,ids].mean(1)

  fig = plt.figure(figsize=(14,6*nrows))
  # fig.suptitle('Stomach_marker_list_V1')
  # Plot cluster
  ax = fig.add_subplot(nrows, 2, 1);                  sc.pl.umap(adataSubGroup, legend_loc='on data', ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format(cluster_key))
  ax = fig.add_subplot(nrows, 2, 2, projection='3d'); sc.pl.umap(adataSubGroup                      , ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format(cluster_key))

  # Plots mean marker genes
  ax = fig.add_subplot(nrows, 2, 3);                  sc.pl.umap(adataSubGroup, legend_loc=None     , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))
  ax = fig.add_subplot(nrows, 2, 4, projection='3d'); sc.pl.umap(adataSubGroup                      , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))

  # Plot individual marker genes
  m=5; n=6
  for i,mgene in enumerate(validgenes):
    # print(i+m, i+n, mgene)
    ax = fig.add_subplot(nrows, 2, i+m);                  sc.pl.umap(adataSubGroup, legend_loc=None     , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows, 2, i+n, projection='3d'); sc.pl.umap(adataSubGroup                      , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
    m+=1; n+=1

  plt.tight_layout()
  plt.savefig("{0}/32_{1}_mouse_cellatlas_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, k) , bbox_inches='tight', dpi=100); plt.close('all')

# 9.5) Save the cellType assigned adataSubGroup into a file
# Write the adataSubGroup and cadataSubGroup object to file
adataSubGroupfile  = "{0}/06_markerGenes_{1}_adataSubGroup.h5ad" .format(dataDir, projName); adataSubGroup.write(adataSubGroupfile)
# # Read back the corrected adataSubGroup object
# adataSubGroupfile  = "{0}/06_markerGenes_{1}_adataSubGroup.h5ad" .format(dataDir, projName); markeradataSubGroup  = sc.read_h5ad(adataSubGroupfile)

################################################################################
# adataSubGroup = markeradataSubGroup.copy()

# 10) subCellType assignments
# Subcluster keys
cluster_key            = "louvain_r0.5"
cluster_bname          = "louvain_r05"
subplot_title_fontsize = 12
subplot_title_width    = 50

# 10.1) Add new categories 
# Categories to rename
adataSubGroup.obs[cluster_key].cat.categories
# Get a new cell type column from the annotation of the louvain_r0.5 clusters
adataSubGroup.obs['subCellType'] = adataSubGroup.obs[cluster_key]
# Add new categories
adataSubGroup.obs['subCellType'].cat.add_categories(['Unknown0','Unknown1','Unknown2','Unknown3','Unknown4','ECL_Restin_Macrophages','ProgenitorAtNeck_PitMucous','Tumor'], inplace=True) 
# Get a new subcluster column
# 0           = 'Unknown0'
# 1           = 'Unknown1'
# 2           = 'Unknown2'
# 3           = 'Unknown3'
# 5           = 'Unknown4'
# 7           = 'ECL_Restin_Macrophages'
# 8           = 'ProgenitorAtNeck_PitMucous'
# 4,6,8,9,10  = 'Tumor'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='0']  = 'Unknown0'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='1']  = 'Unknown1'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='2']  = 'Unknown2'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='3']  = 'Unknown3'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='5']  = 'Unknown4'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='7']  = 'ECL_Restin_Macrophages'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='8']  = 'ProgenitorAtNeck_PitMucous'
adataSubGroup.obs['subCellType'].loc[(adataSubGroup.obs['subCellType']=='4')|(adataSubGroup.obs['subCellType']=='5')|(adataSubGroup.obs['subCellType']=='6')|(adataSubGroup.obs['subCellType']=='7')|(adataSubGroup.obs['subCellType']=='9')|(adataSubGroup.obs['subCellType']=='10')]  = 'Tumor'
# Remove old categories
adataSubGroup.obs['subCellType'].cat.remove_categories(adataSubGroup.obs[cluster_key].cat.categories.tolist(), inplace=True)
# List new categories
adataSubGroup.obs['subCellType'].cat.categories
# Draw Umaps with the categories
# Subcluster keys
cluster_key            = "subCellType"
cluster_bname          = "subCellType"
subplot_title_fontsize = 12
subplot_title_width    = 50

# Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
cluster_key_groups = adataSubGroup.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = adataSubGroup.obs[cluster_key].value_counts().to_dict()

def adjust_title(ax):
  '''https://stackoverflow.com/questions/55197674/matplotlib-prevent-subplot-title-from-being-wider-than-subplot'''
  title = ax.title
  ax.figure.canvas.draw()
  def _get_t():
      ax_width = ax.get_window_extent().width
      ti_width = title.get_window_extent().width
      return ax_width/ti_width

  while _get_t() <= 1 and title.get_fontsize() > 1:        
      title.set_fontsize(title.get_fontsize()-1)


# Louvain UMAPs
ncols  = len(cluster_key_groups) + 1
fig = plt.figure(figsize=(8*ncols, 14))
fig.suptitle("{0} UMAP".format(cluster_key))
# Main Louvain Cluster
ax = fig.add_subplot(2, ncols, 1); sc.pl.umap(adataSubGroup, legend_loc='on data', ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, show=False); adjust_title(ax)
ax = fig.add_subplot(2, ncols, 2, projection='3d'); sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False); adjust_title(ax)
# Partial visualizaton of a subset of groups in embedding
m=3; n=4
for i,b in enumerate(cluster_key_groups):
  print(i, b)
  ax = fig.add_subplot(2, ncols, i+m);                  sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b])); adjust_title(ax)
  ax = fig.add_subplot(2, ncols, i+n, projection='3d'); sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9,  projection='3d', show=False); ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b])); adjust_title(ax)
  m+=1; n+=1

plt.tight_layout()
plt.savefig("{0}/06_subGroup_{1}_clustering_subCellType_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=200); plt.close('all')

# 10.2) Plot separate bar plots, coloured in by cluster annotation, for each tissue
# Convert palette into colormap
clcmap = ListedColormap(adataSubGroup.uns['subCellType_colors'])
# Get the DF of tissue and clusters
clusterBatchDF = adataSubGroup.obs[['batch','subCellType']].copy()
# Replace batch number with batch names
clusterBatchDF.replace({'batch': tissueIdDict}, inplace=True)
# Remove index for groupby
clusterBatchDF.reset_index(drop=True, inplace=True)
# Remove unused categories (0,1,..)
clusterBatchDF['subCellType'] = clusterBatchDF.subCellType.cat.remove_unused_categories()
# Get the number of cells for each cluster in every tissue
ncellsClusterBatchDF = clusterBatchDF.groupby(['batch','subCellType']).size()
# Get the percent of cells for each cluster in every tissue 
pcellsClusterBatchDF = pd.crosstab(index=clusterBatchDF['batch'], columns=clusterBatchDF['subCellType'], values=clusterBatchDF['subCellType'], aggfunc='count', normalize='index')

# Plot the barplots
fig = plt.figure(figsize=(16,6)); fig.suptitle("Cells for each {0} in each tissue".format('subCellType'))
# plot numbers of cells
ax = fig.add_subplot(1, 2, 1); ncellsClusterBatchDF.unstack().plot(kind='barh', stacked=True, colormap=clcmap, ax=ax, legend=None, title="Number of cells", edgecolor=None)
# plot percent of cells
ax = fig.add_subplot(1, 2, 2); pcellsClusterBatchDF.plot(kind='barh',stacked=True, colormap=clcmap, ax=ax, title="% of cells")
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='subCellType', title_fontsize=12)
plt.tight_layout() # For non-overlaping subplots
plt.savefig("{0}/06_subGroup_{1}_clustering_{2}_tissueID_subCellType_barplot.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 10.3) Calculate marker genes (one vs. rest)
sc.tl.rank_genes_groups(adataSubGroup, groupby='subCellType', key_added='rank_genes')
# Plot marker genes
sc.pl.rank_genes_groups(adataSubGroup, key='rank_genes', fontsize=12, show=False)
plt.savefig("{0}/06_{1}_{2}_marker_genes_ranking_subCellType.png".format(plotsDir, bname, 'subCellType') , bbox_inches='tight', dpi=175); plt.close('all')

# 10.4) Calculate pairwise marker genes list
# Get all subCellTypes into the list
subCellTypeCategories = adataSubGroup.obs['subCellType'].cat.categories.tolist()
# Get the list of all unique pairwise combinations
subCellTypePairComb = [comb for comb in combinations(subCellTypeCategories, 2)]
# Get pairwise plots dir
pwsgPlotsDir = "{0}/subGroup/rankedGenes".format(plotsDir); create_dir(pwsgPlotsDir)
# Get pairwise data dir
pwsgDataDir  = "{0}/subGroup/rankedGenes".format(dataDir); create_dir(pwsgDataDir)
# Calculate pairwise marker genes  
for grp,ref in subCellTypePairComb:
  print("- Calculating pairwise marker genes for group_v_reference: {0}_v_{1}".format(grp, ref))
  # Get genes ranking
  keyName = 'rank_genes_{0}_v_{1}'.format(grp,ref)
  sc.tl.rank_genes_groups(adataSubGroup, groupby='subCellType', groups=[grp], key_added=keyName, reference=ref, n_genes=adataSubGroup.shape[1])
  # Plot top 20 ranked genes
  sc.pl.rank_genes_groups(adataSubGroup, key=keyName, groups=[grp], fontsize=12, show=False)
  # Save it in a figure
  plt.savefig("{0}/06_{1}_all_subCellType_{2}_v_{3}.png".format(pwsgPlotsDir, bname, grp, ref) , bbox_inches='tight'); plt.close('all')
  # Get the dataframe of DE parameters
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(adataSubGroup.uns[keyName][n])[grp]
  # Add treatment and reference group name
  ngDF['Treatment'] = grp
  ngDF['Reference'] = ref
  # Save the dataframe
  ngDF.to_csv("{0}/04_{1}_{2}.txt".format(pwsgDataDir, projName, keyName), sep='\t', header=True, index=False, float_format='%.2g')

# 10.5) Save the coordinates of the umap
rowIdx     = adataSubGroup.obs.index.values.tolist()
umapCordDF = pd.DataFrame(adataSubGroup.obsm['X_umap'],index=rowIdx, columns=['UMAP1','UMAP2', 'UMAP3'])
umapCordDF.to_csv("{0}/07_{1}_subCellType_assigned_{2}_UMAP_Coordinates.txt" .format(dataDir, projName, cluster_key), sep='\t', header=True, index=True, index_label="CellID")

# 10.7) Save the subCellType assigned adataSubGroup into a file
# Write the adataSubGroup and cadataSubGroup object to file
adataSubGroupfile  = "{0}/07_subCellType_assigned_{1}_adataSubGroup.h5ad" .format(dataDir, projName); adataSubGroup.write(adataSubGroupfile)

############################################################
# 11) Analysis of Tumor cluster in SubGroup
############################################################
subCellTypeadata = adataSubGroup.copy() # (5885, 11023)
# adataSubGroup = subCellTypeadata.copy() # (5885, 11023)

# 8.1) Subcluster Tumors
sc.tl.louvain(adataSubGroup, restrict_to=('subCellType', ['Tumor']), resolution=0.3, key_added='subCellType_Tumor_sub')
#Show the new clustering
if 'subCellType_Tumor_sub_colors' in adataSubGroup.uns:
    del adataSubGroup.uns['subCellType_Tumor_sub_colors']

# Subcluster keys
cluster_key            = "subCellType_Tumor_sub"
cluster_bname          = "subCellType_Tumor_sub"
subplot_title_fontsize = 12
subplot_title_width    = 50

# Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
cluster_key_groups = adataSubGroup.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = adataSubGroup.obs[cluster_key].value_counts().to_dict()

# Louvain UMAPs
ncols  = len(cluster_key_groups) + 1
fig = plt.figure(figsize=(7*ncols, 14))
fig.suptitle("{0} UMAP".format(cluster_key))
# Main Louvain Cluster
ax = fig.add_subplot(2, ncols, 1);                  sc.pl.umap(adataSubGroup, legend_loc='on data', ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(2, ncols, 2, projection='3d'); sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
# Partial visualizaton of a subset of groups in embedding
m=3; n=4
for i,b in enumerate(cluster_key_groups):
  print(i, b)
  ax = fig.add_subplot(2, ncols, i+m);                  sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
  ax = fig.add_subplot(2, ncols, i+n, projection='3d'); sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
  m+=1; n+=1

plt.tight_layout()
plt.savefig("{0}/05_subGroup_{1}_clustering_{2}_UMAP.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Get the new marker genes (one. vs. rest)
sc.tl.rank_genes_groups(adataSubGroup, groupby='subCellType_Tumor_sub', key_added='rank_genes_subCellType_Tumor_sub')
# Plot the new marker genes
sc.pl.rank_genes_groups(adataSubGroup, key='rank_genes_subCellType_Tumor_sub', groups=['Tumor,0','Tumor,1','Tumor,2', 'Tumor,3'], fontsize=12, show=False)
plt.savefig("{0}/06_subGroup_{1}_clustering_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 10.4) Calculate pairwise marker genes list
# Get all subCellTypes into the list
subCellTypeCategories = adataSubGroup.obs[cluster_key].cat.categories.tolist()
# Get the list of all unique pairwise combinations
subCellTypePairComb = [comb for comb in combinations(subCellTypeCategories, 2)]
# Get pairwise plots dir
tspwsgPlotsDir = "{0}/subGroup/rankedGenes/tumorSubCluster".format(plotsDir); create_dir(tspwsgPlotsDir)
# Get pairwise data dir
tspwsgDataDir  = "{0}/subGroup/rankedGenes/tumorSubCluster".format(dataDir); create_dir(tspwsgDataDir)
# Calculate pairwise marker genes  
for grp,ref in subCellTypePairComb:
  print("- Calculating pairwise marker genes for group_v_reference: {0}_v_{1}".format(grp, ref))
  # Get genes ranking
  keyName = 'rank_genes_{0}_v_{1}'.format(grp,ref)
  sc.tl.rank_genes_groups(adataSubGroup, groupby=cluster_key, groups=[grp], key_added=keyName, reference=ref, n_genes=adataSubGroup.shape[1])
  # Plot top 20 ranked genes
  sc.pl.rank_genes_groups(adataSubGroup, key=keyName, groups=[grp], fontsize=12, show=False)
  # Save it in a figure
  plt.savefig("{0}/04_{1}_all_subCellType_{2}_v_{3}.png".format(tspwsgPlotsDir, bname, grp, ref) , bbox_inches='tight'); plt.close('all')
  # Get the dataframe of DE parameters
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(adataSubGroup.uns[keyName][n])[grp]
  # Add treatment and reference group name
  ngDF['Treatment'] = grp
  ngDF['Reference'] = ref
  # Save the dataframe
  ngDF.to_csv("{0}/04_{1}_{2}.txt".format(tspwsgDataDir, projName, keyName), sep='\t', header=True, index=False, float_format='%.2g')

# 10.5) Save the normalized, log transformed, batch and cell cycle corrected data
CellTypeDF                                = adata.to_df()
CellTypeDF['TissueID']                    = adata.obs['tissueID']
CellTypeDF['OriginalLouvainCluster']      = adata.obs['louvain_r0.5']
CellTypeDF['OriginalCellType']            = adata.obs['cellType']
CellTypeDF['SubCellLouvainCluster']       = adataSubGroup.obs['louvain_r0.5']
CellTypeDF['SubCellType']                 = adataSubGroup.obs['subCellType']
CellTypeDF['SubCellTypeTumorSubClusters'] = adataSubGroup.obs['subCellType_Tumor_sub']
CellTypeDF.to_csv("{0}/07_{1}_scran_normalized_counts_annotation.txt" .format(dataDir, projName, cluster_key), sep='\t', header=True, index=True, index_label="CellID")

# subCellTypeDF.T.to_csv("{0}/04_normalizedRaw_{1}_customClusters.txt".format(countsDir, projName), sep='\t', header=True, index=True, index_label="GeneSymbol", float_format='%.2g')

# 10.6) Save the subCellType assigned adataSubGroup into a file
# Write the adataSubGroup and cadataSubGroup object to file
adataSubGroupfile  = "{0}/08_{1}_subCellType_tumorSubCluster_adataSubGroup.h5ad" .format(dataDir, projName); adataSubGroup.write(adataSubGroupfile)

# Finished on 2020-04Apr-22 19:55
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

