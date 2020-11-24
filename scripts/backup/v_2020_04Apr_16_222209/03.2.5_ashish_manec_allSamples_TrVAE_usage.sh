# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

# Source: https://nbviewer.jupyter.org/github/theislab/trVAE/blob/master/examples/trVAE_Haber.ipynb

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

# For X11 display
import matplotlib
# matplotlib.use('TkAgg')
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt


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
minGenesPerCell = 25
minCountPerCell = 150
maxCountPerCell = 50000 
minCellsPergene = 25
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
# Memory usage: current 0.89 GB, difference +0.89 GB
# Running Scanpy 1.4.5.1, on 2020-04-15 13:15.
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
                    '0':'bulk1001', 
                    '1':'bulk997', 
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
minGenesPerCell = 25
minCountPerCell = 150
maxCountPerCell = 50000 
minCellsPergene = 25

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
# Number of cells after gene filter: 5885

# 2.6) Filter genes according to identified QC thresholds:
# Min minCellsPergene cells - filters out 0 count genes
print('Total number of genes: {:d}'.format(adata.n_vars))
sc.pp.filter_genes(adata, min_cells=minCellsPergene)
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
# Total number of genes: 11029
# Number of genes after cell filter: 11023

# 2.7) Compute highly variable genes (HVGs)
# Next, we first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000, batch_key = 'tissueID')
print("Highly variable genes intersection: %d"%sum(adata.var.highly_variable_intersection))
print("Number of batches where gene is variable:")
print(adata.var.highly_variable_nbatches.value_counts())
# Highly variable genes intersection: 983
# Number of batches where gene is variable:
# 1    3460
# 0    2978
# 2    2200
# 3    1402
# 4     983
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

#########################################################################
normadata = adata.copy()
rawadatas = adatas.copy()
# 5) Technical correction: Batch Correction using TrVAE

# 5.1) Normalizing & Extracting Top 1000 Highly Variable Genes
# We can preserve more genes (up to 7000 like scGen) but in order to train the network quickly, we will extract top 1000 genes
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=1000)
adata = adata[:, adata.var['highly_variable']]
adata

# 5.2) Train/Test Split
train_adata, valid_adata = trvae.utils.train_test_split(adata, train_frac=0.80)
train_adata.shape, valid_adata.shape
# ((4708, 1000), (1177, 1000))

# 5.3) Calculate number of batches
n_conditions = len(train_adata.obs[condition_key].unique().tolist())
conditions = adata.obs[condition_key].unique().tolist()
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
# NOTE: This will take 8 minutes on WS4 with single CPU
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
              )

# 5.6) Getting batch-corrected adata
# Now two matrices have been added to adata
# mmd_latent: (numpy ndarray) output of MMD Layer in trVAE
# reconstructed: (numpy ndarray) reconstructed data with dimension of original feature space
# z_latent: (numpy ndarray) output of bottleneck layer of trVAE (optional)
# For evaluating how good trVAE has corrected the batches, we recommend using mmd_latent matrix.
labels, _ = trvae.tl.label_encoder(adata, condition_key=condition_key, label_encoder=condition_encoder)
network.get_corrected(adata, labels, return_z=True)

# 5.7) MMD Layer UMAP visualization
mmd_latent = adata.obsm['mmd_latent']
mmd_adata = sc.AnnData(mmd_latent, obs=adata.obs)
mmd_adata
sc.pp.neighbors(adata, random_state = 2105, n_neighbors=10, use_rep = "mmd_latent")
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

############################################################
############################################################
# Using TrVAE for batch correction
############################################################
trvaeBCadata = adata.copy() # (5885, 11023)

# 6) Plot UMAP for all the tumor cells 
# Color the cells that have human myc and ires
cellBarCodes = pd.read_csv('/media/rad/HDD2/temp_manec/hgMycIresCd2_cellIDs.txt', sep="\t", header=None).values.tolist()
cellBarCodes = pd.read_csv('/media/rad/HDD2/temp_manec/hgMycIresCd2_humanCd2_cellIDs.txt', sep="\t", header=None).values.tolist()
cl  = sum(cellBarCodes, [])
ucl = get_unique_list(cl)
# In [34]: len(ucl)
# Out[34]: 1743

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
# 1     692
# 2     685
# 3     653
# 4     650
# 5     571
# 6     541
# 7     475
# 8     434
# 9     401
# 10    370
# 11    325
# 12    218
# 13    203
# 14    157
# 15    137
# 16    130
# 17    114

# 4.3) Visualizations

# Calculations for the visualizations
# sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
# sc.pp.neighbors(adata, random_state = 2105, use_rep = "SC")
# sc.tl.umap(adata, random_state = 2105, n_components=3)

# Plot visualizations
# Visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=["{0}".format(cluster_key), 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_normTrVAE_{1}_clustering_all_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adata, color=["{0}".format(cluster_key), 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
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

# 7.2) Marker genes & cluster annotation
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby=cluster_key, key_added='rank_genes')

# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes', fontsize=12, show=False)
plt.savefig("{0}/03_normTrVAE_{1}_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# Annotation of cluster r_0.5 with known marker genes
markerDir = "{0}/markerDir".format(plotsDir); create_dir(markerDir)

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

# Generate the UMAPs for each marker categories
for k in marker_genes.keys():
  ids = np.in1d(adata.var_names, marker_genes[k])
  adata.obs['{0}_marker_expr'.format(k)] = adata.X[:,ids].mean(1)
  sc.pl.umap(adata, color=["{0}".format(cluster_key),'{0}_marker_expr'.format(k)], color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9,  projection='3d', show=False)
  plt.savefig("{0}/31_{1}_marker_genes_stomach_{2}_UMAPs_3D.png".format(markerDir, bname, k) , bbox_inches='tight', dpi=175); plt.close('all')

# Generate the UMAPs for each marker categories
plt.figure(figsize=(100,8))
for k in ma_marker_genes.keys():
  ids = np.in1d(adata.var_names, ma_marker_genes[k])
  adata.obs['{0}_ma_marker_expr'.format(k)] = adata.X[:,ids].mean(1)
  sc.pl.umap(adata, color=["{0}".format(cluster_key),'{0}_ma_marker_expr'.format(k)], color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9,  projection='3d', show=False)
  plt.savefig("{0}/32_{1}_mouse_cellatlas_marker_genes_stomach_{2}_UMAPs_3D.png".format(markerDir, bname, k) , bbox_inches='tight', dpi=175); plt.close('all')

# Plot Final Marker genes
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby=cluster_key, key_added='rank_genes_louvain')
# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_louvain', fontsize=12, show=False)
plt.savefig("{0}/{1}_louvain_marker_genes_ranking.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Erythrocytes
sc.pl.umap(adata, color=["{0}".format(cluster_key),'Hbb-bs', 'Hba-a1','Hba-a2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Erythrocytes') , bbox_inches='tight', dpi=175); plt.close('all')

# Restin
sc.pl.umap(adata, color=["{0}".format(cluster_key),'S100a8','S100a9'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Restin') , bbox_inches='tight', dpi=175); plt.close('all')

# Tcells
sc.pl.umap(adata, color=["{0}".format(cluster_key),'Cd3d'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Tcells') , bbox_inches='tight', dpi=175); plt.close('all')

# Endothelial/Epithelial_Igfbp3⁺
sc.pl.umap(adata, color=["{0}".format(cluster_key),'Egfl7','Sparc', 'Col4a1','Plvap', 'Ifitm3','Esam', 'Cdh5', 'Igfbp3','Plpp3','Kdr'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Endothelial_Epithelial_Igfbp3pos') , bbox_inches='tight', dpi=175); plt.close('all')

# Endothelial
sc.pl.umap(adata, color=["{0}".format(cluster_key),'Plvap', 'Ctla2a','Eng'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Endothelial') , bbox_inches='tight', dpi=175); plt.close('all')

# # Pancreas (Identified from louvain_r15_marker_genes_ranking)
# # This particular plot is providing no useful information
# # Check: https://www.proteinatlas.org/ENSG00000125691-RPL23/summary/rna
# sc.pl.umap(adata, color=["{0}".format(cluster_key),'Tpt1','Eef1b2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
# plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Pancreas') , bbox_inches='tight', dpi=175); plt.close('all')

# Save the louvain information in external file
louvainsDF = pd.DataFrame(adata.obs[cluster_key])
louvainsDF.to_csv("{0}/03_{1}_louvains.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="cellId")

# Finished on 2020-04Apr-15 17:51 am
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#---------------------------------------------------------------------
# Categories to rename
adata.obs[cluster_key].cat.categories

# Get a new cell type column from the annotation of the louvain_r0.5 clusters
adata.obs['cellType'] = adata.obs[cluster_key]

# Add new categories
adata.obs['cellType'].cat.add_categories(['Birc5⁺/basal', 'Tcells', 'Dendritic', 'Macrophages','Endothelial', 'Endothelial/Epithelial_Igfbp3⁺', 'Erythrocytes', 'Fibroblasts','Restin_like_gamma', 'Tumor', 'Pit_cells', 'Parietal', 'Pancreas'], inplace=True) 

# Get a new subcluster column
# 0           = 'Birc5⁺/basal'
# 1           = 'Tcells'
# 2           = 'Dendritic'
# 3,7,8,11,12 = 'Tumor'
# 4,6         = 'Pit_cells'
# 5           = 'Erythrocytes'
# 9           = 'Macrophages'
# 10          = 'Endothelial'
# 13          = 'Fibroblasts'
# 14          = 'Pancreas'
# 15          = 'Restin_like_gamma'
# 16          = 'Endothelial/Epithelial_Igfbp3⁺'
# 17          = 'Parietal'

adata.obs['cellType'].loc[adata.obs['cellType' ]=='0' ]  = 'Birc5⁺/basal'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='1' ]  = 'Tcells'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='2' ]  = 'Dendritic'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='5' ]  = 'Erythrocytes'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='9' ]  = 'Macrophages'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='10']  = 'Endothelial'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='13']  = 'Fibroblasts'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='14']  = 'Pancreas'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='15']  = 'Restin_like_gamma'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='16']  = 'Endothelial/Epithelial_Igfbp3⁺'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='17']  = 'Parietal'
adata.obs['cellType'].loc[(adata.obs['cellType']=='6')|(adata.obs['cellType']=='4')]  = 'Pit_cells'
adata.obs['cellType'].loc[(adata.obs['cellType']=='3')|(adata.obs['cellType']=='7')|(adata.obs['cellType']=='8')|(adata.obs['cellType']=='11')|(adata.obs['cellType']=='12')]  = 'Tumor'

# Remove old categories
adata.obs['cellType'].cat.remove_categories(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'], inplace=True)

# List new categories
adata.obs['cellType'].cat.categories

# Draw Umaps with the categories
sc.pl.umap(adata, color=['cellType'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_normScanorama_{1}_clustering_cellType_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.umap(adata, color=['cellType'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, legend_loc='on data', show=False)
plt.savefig("{0}/03_normScanorama_{1}_clustering_cellType_legend_onData_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

plt.figure(figsize=(50,50))
sc.pl.umap(adata, color=['cellType'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, legend_loc='on data', show=False,zorder=0)
plt.savefig("{0}/03_normScanorama_{1}_clustering_cellType_legend_onData_UMAP.jpg".format(plotsDir, bname) , bbox_inches='tight', dpi=600, rasterized=True, edgecolor='white'); plt.close('all')
