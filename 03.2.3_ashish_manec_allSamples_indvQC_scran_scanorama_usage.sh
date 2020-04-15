# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

# NOTE: plt.subplot(number_of_rows, number_of_columns, index_of_the_subplot) 
ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)

# Loading the python libraries
import scanpy as sc
import scanpy.external as sce
from gjainPyLib import *
import pickle
import logging
import scanorama

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
projName        = "alltissues_except1079" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/{0}".format(projName); create_dir("{0}".format(output_dir))
ccGenes_macosko = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/macosko_cell_cycle_genes_mmu.txt"
ccGenes_regev   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/regev_lab_cell_cycle_genes_mmu.txt"
minGenesPerCell = 50
minCountPerCell = 100
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
# Running Scanpy 1.4.5.1, on 2020-04-14 22:53.
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
minCountPerCell = 25
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

# 1.6) Checking the total size of the data set
adata.shape # We have 5887 cells and 6068 genes in the dataset

# Reset the parameters for merged data
minGenesPerCell = 50
minCountPerCell = 100
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
# Min 5 cells - filters out 0 count genes
print('Total number of genes: {:d}'.format(adata.n_vars))
sc.pp.filter_genes(adata, min_cells=minCellsPergene)
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
# Number of genes after cell filter: 11023
# Highly variable genes intersection: 983

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

# 5) Technical correction: Batch Correction using Scanorama
# 5.1) Detect variable genes
# # As the stored AnnData object contains scaled data based on variable genes, we need to make a new object with the raw counts and normalized it again. Variable gene selection should not be performed on the scaled data object, only do normalization and log transformation before variable genes selection.
# adata2 = sc.AnnData(X=adata.raw.X, var=adata.raw.var, obs = adata.obs)
# sc.pp.normalize_per_cell(adata2, counts_per_cell_after=1e4)
# sc.pp.log1p(adata2)

adata2 = normadata.copy()

# Detect variable genes for the full dataset
sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Highly variable genes: %d"%sum(adata2.var.highly_variable))
var_genes_all = adata2.var.highly_variable
# Highly variable genes: 2384

# Detect variable genes in each dataset separately using the batch_key parameter.
sc.pp.highly_variable_genes(adata2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'tissueID')
var_genes_batch = adata2.var.highly_variable_nbatches > 0
print("Any batch var genes: %d"%sum(var_genes_batch))
print("All data var genes: %d"%sum(var_genes_all))
print("Overlap: %d"%sum(var_genes_batch & var_genes_all))
print("Variable genes in all batches: %d"%sum(adata2.var.highly_variable_nbatches ==3))
print("Overlap batch instersection and all: %d"%sum(var_genes_all & adata2.var.highly_variable_intersection))

# Select all genes that are variable in at least 2 datasets and use for remaining analysis.
var_select = adata2.var.highly_variable_nbatches > 2
var_genes = var_select.index[var_select]

# 5.2 )Data integration
# Split per batch into new objects.
batches = ['bulk997','bulk1001','bulk1018','stomach1001']
alldata = {}
for batch in batches:
    alldata[batch] = adata2[adata2.obs['tissueID'] == batch,]

# Subset the individual dataset to the variable genes
alldata2 = dict()
for ds in alldata.keys():
    print(ds)
    alldata2[ds] = alldata[ds][:,var_genes]

# Convert to list of AnnData objects
adatas = list(alldata2.values())

# Run scanorama.integrate
scanorama  = scanorama.integrate_scanpy(adatas, dimred = 50,)

# Returns a list of 4 np.ndarrays with 50 columns.
print(scanorama[0].shape)
print(scanorama[1].shape)
print(scanorama[2].shape)
print(scanorama[3].shape)
# (971, 50)
# (3462, 50)
# (671, 50)
# (781, 50)

# Make into one matrix.
all_s = np.concatenate(scanorama)
print(all_s.shape) # (5885, 50)

# Add to the AnnData object
adata.obsm["SC"] = all_s

# Calculations for the visualizations
sc.pp.neighbors(adata, random_state = 2105, use_rep = "SC")
sc.tl.umap(adata, random_state = 2105, n_components=3)

# Plot the UMAPs
fig = plt.figure(figsize=(16,13))
# 2D projection
ax = fig.add_subplot(2, 2, 1);                  sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 2, 2);                  sc.pl.umap(adata, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Scran Scanorama UMAP")
# 3D projection
ax = fig.add_subplot(2, 2, 3, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 2, 4, projection='3d'); sc.pl.umap(adata                , ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Scran Scanorama UMAP")
plt.tight_layout()
plt.savefig("{0}/03_norm_all_batchCorrection_{1}_tissueID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

############################################################
############################################################
# Using SCANORAMA for batch correction
############################################################
scanoramaBCadata = adata.copy() # (5885, 11023)

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
sc.pp.neighbors(adata, random_state = 2105, use_rep = "SC")
sc.tl.umap(adata, random_state = 2105, n_components=3)

# Plot visualizations
# Visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_normScanorama_{1}_clustering_all_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/03_normScanorama_{1}_clustering_all_louvain_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

fig = plt.figure(figsize=(24,8))
# 2D projection
ax = fig.add_subplot(2, 4, 1);                  sc.pl.umap(adata             ,                  ax=ax, color="louvain_r0.5", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="louvain_r0.5 UMAP")
ax = fig.add_subplot(2, 4, 2);                  sc.pl.umap(adata   , ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="tissueID UMAP")
ax = fig.add_subplot(2, 4, 3);                  sc.pl.umap(adata         , legend_loc=None, ax=ax, color="log_counts"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="log_counts UMAP")
ax = fig.add_subplot(2, 4, 4);                  sc.pl.umap(adata, legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="mt_frac UMAP")
# 3D projection
ax = fig.add_subplot(2, 4, 5, projection='3d'); sc.pl.umap(adata             , legend_loc=None, ax=ax, color="louvain_r0.5", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="louvain_r0.5 UMAP")
ax = fig.add_subplot(2, 4, 6, projection='3d'); sc.pl.umap(adata   , ax=ax, color="tissueID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="tissueID UMAP")
ax = fig.add_subplot(2, 4, 7, projection='3d'); sc.pl.umap(adata         , legend_loc=None, ax=ax, color="log_counts"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="log_counts UMAP")
ax = fig.add_subplot(2, 4, 8, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="mt_frac UMAP")
plt.tight_layout()
plt.savefig("{0}/03_normScanorama_{1}_louvain_tissueID_counts_mtfrac_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Louvain UMAPs
fig = plt.figure(figsize=(16,6))
fig.suptitle('louvain')
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  
sc.pl.umap(adata, legend_loc=None, ax=ax, color="louvain", palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(1, 2, 2, projection='3d'); 
sc.pl.umap(adata, ax=ax, color="louvain", palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/03_normScanorama_{1}_clustering_louvain_r05_UMAP_2D3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 7.2) Marker genes & cluster annotation
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby='louvain', key_added='rank_genes')

# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes', fontsize=12, show=False)
plt.savefig("{0}/03_normScanorama_{1}_louvain_r05_marker_genes_ranking.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

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
  sc.pl.umap(adata, color=['louvain','{0}_marker_expr'.format(k)], color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9,  projection='3d', show=False)
  plt.savefig("{0}/31_{1}_marker_genes_stomach_{2}_UMAPs_3D.png".format(markerDir, bname, k) , bbox_inches='tight', dpi=175); plt.close('all')

# Generate the UMAPs for each marker categories
plt.figure(figsize=(100,8))
for k in ma_marker_genes.keys():
  ids = np.in1d(adata.var_names, ma_marker_genes[k])
  adata.obs['{0}_ma_marker_expr'.format(k)] = adata.X[:,ids].mean(1)
  sc.pl.umap(adata, color=['louvain','{0}_ma_marker_expr'.format(k)], color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9,  projection='3d', show=False)
  plt.savefig("{0}/32_{1}_mouse_cellatlas_marker_genes_stomach_{2}_UMAPs_3D.png".format(markerDir, bname, k) , bbox_inches='tight', dpi=175); plt.close('all')

# Plot Final Marker genes
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby='louvain', key_added='rank_genes_louvain')
# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_louvain', fontsize=12, show=False)
plt.savefig("{0}/{1}_louvain_marker_genes_ranking.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Erythrocytes
sc.pl.umap(adata, color=['louvain','Hbb-bs', 'Hba-a1','Hba-a2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Erythrocytes') , bbox_inches='tight', dpi=175); plt.close('all')

# Restin
sc.pl.umap(adata, color=['louvain','S100a8','S100a9'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Restin') , bbox_inches='tight', dpi=175); plt.close('all')

# Tcells
sc.pl.umap(adata, color=['louvain','Cd3d'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Tcells') , bbox_inches='tight', dpi=175); plt.close('all')

# Endothelial/Epithelial_Igfbp3⁺
sc.pl.umap(adata, color=['louvain','Egfl7','Sparc', 'Col4a1','Plvap', 'Cd93','Ifitm3','Esam', 'Cdh5', 'Igfbp3','Plpp3','Kdr','Sptbn1'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Endothelial_Epithelial_Igfbp3pos') , bbox_inches='tight', dpi=175); plt.close('all')

# Endothelial
sc.pl.umap(adata, color=['louvain','Plvap', 'Cd34', 'Ctla2a','Cd93','Ramp2','Eng'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Endothelial') , bbox_inches='tight', dpi=175); plt.close('all')

# Pancreas (Identified from louvain_r15_marker_genes_ranking)
# This particular plot is providing no useful information
# Check: https://www.proteinatlas.org/ENSG00000125691-RPL23/summary/rna
sc.pl.umap(adata, color=['louvain','Rps24','Rps11','Rps3','Rpl23','Tpt1','Eef1b2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Pancreas') , bbox_inches='tight', dpi=175); plt.close('all')

# Save the louvain information in external file
louvainsDF = pd.DataFrame(adata.obs['louvain'])
louvainsDF.to_csv("{0}/03_{1}_louvains.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="cellId")





# Finished on 2020-04Apr-14
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# --------------------------------------------------------
# Plot UMAP for all the tumor cells 
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
# sc.pl.umap(adata, color='hgMycIresCd2', use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
# plt.savefig("{0}/03_normScanorama_{1}_Tumor_hgMycIresCd2_CellIDs_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
fig = plt.figure(figsize=(16,6))
fig.suptitle('hgMycIresCd2')
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  
sc.pl.umap(adata, legend_loc=None, ax=ax, color="hgMycIresCd2", color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(1, 2, 2, projection='3d'); 
sc.pl.umap(adata, ax=ax, color="hgMycIresCd2", color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/03_normScanorama_{1}_Tumor_hgMycIresCd2_CellIDs_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# # --------------------------------------------------------
# # Get tumors for all individual vectors ('humanMyc', 'gap', 'ires', 'humanCd2')
# # Unique cell barcode list of list
# ucblofl = list()
# ucbd    = defaultdict(list)
# for t in ['humanMyc', 'gap', 'ires', 'humanCd2']:
#   # Cell barcode data frame
#   cbDF = pd.read_csv('/media/rad/HDD2/temp_manec/hgMycIresCd2_{0}_cellIDs.txt'.format(t), sep="\t", header=None).values.tolist()

#   # Unique cell barcode list
#   ucbl = get_unique_list(sum(cbDF, []))
#   ucbd[t] = ucbl
#   ucblofl.append(ucbl)

# import upsetplot
# from upsetplot import from_memberships, from_contents
# upsetContent = from_contents(ucbd)
# memberships  = from_memberships(ucblofl)

# fig = plt.figure(figsize=(15,5), dpi=175)
# upsetplot.plot(upsetContent, sum_over=False,show_counts='%d', fig=fig)
# plt.savefig("{0}/03_normScanorama_{1}_Tumor_hgMycIresCd2_individual_CellIDs_upsetPlot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# # Plot only selected groups
# upcDF        = upsetContent.query('humanMyc == True and gap == False and ires == False and humanCd2 == False')
# cellBarCodes = upcDF.values.tolist()
# upcDF        = upsetContent.query('humanMyc == False and gap == True and ires == False and humanCd2 == False')
# cellBarCodes = upcDF.values.tolist()
# upcDF        = upsetContent.query('humanMyc == False and gap == False and ires == True and humanCd2 == False')
# cellBarCodes = upcDF.values.tolist()
# upcDF        = upsetContent.query('humanMyc == False and gap == False and ires == False and humanCd2 == True')
# cellBarCodes = upcDF.values.tolist()

# cl  = sum(cellBarCodes, [])
# ucl = get_unique_list(cl)
# # In [34]: len(ucl)
# # Out[34]: 1743

# mylist = adata.obs.index.values
# humaniresmyc = list()
# for e in mylist: 
#   flag = 0
#   for s in ucl: 
#       if s in e: 
#           flag = 1 
#           break
#   humaniresmyc.append(flag)

# adata.obs['hgMycIresCd2'] = humaniresmyc
# sc.pl.umap(adata, color='hgMycIresCd2', use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

#---------------------------------------------------------------------
# Categories to rename
adata.obs['louvain_r0.5'].cat.categories

# Get a new cell type column from the annotation of the louvain_r0.5 clusters
adata.obs['cellType'] = adata.obs['louvain_r0.5']

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

# Plot Final Marker genes
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby='cellType', key_added='rank_genes_cellType')
# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_cellType', fontsize=12, show=False)
plt.savefig("{0}/{1}_cellType_marker_genes_ranking.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Erythrocytes
sc.pl.umap(adata, color=['cellType','Hbb-bs', 'Hba-a1','Hba-a2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Erythrocytes') , bbox_inches='tight', dpi=175); plt.close('all')

# Restin
sc.pl.umap(adata, color=['cellType','Retnlg','S100a8','S100a9'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Restin') , bbox_inches='tight', dpi=175); plt.close('all')

# Tcells
sc.pl.umap(adata, color=['cellType','Sh2d1a','Cd3d','Cd3e','Cd8a'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Tcells') , bbox_inches='tight', dpi=175); plt.close('all')

# Endothelial/Epithelial_Igfbp3⁺
sc.pl.umap(adata, color=['cellType','Egfl7','Sparc', 'Col4a1','Plvap', 'Cd93','Ifitm3','Esam', 'Cdh5', 'Igfbp3','Plpp3','Kdr','Sptbn1'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Endothelial_Epithelial_Igfbp3pos') , bbox_inches='tight', dpi=175); plt.close('all')

# Endothelial
sc.pl.umap(adata, color=['cellType','Plvap', 'Cd34', 'Ctla2a','Cd93','Ramp2','Eng'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Endothelial') , bbox_inches='tight', dpi=175); plt.close('all')

# Pancreas (Identified from louvain_r15_marker_genes_ranking)
# This particular plot is providing no useful information
# Check: https://www.proteinatlas.org/ENSG00000125691-RPL23/summary/rna
sc.pl.umap(adata, color=['cellType','Rps24','Rps11','Rps3','Rpl23','Tpt1','Eef1b2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/30_{1}_manually_annotated_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, 'Pancreas') , bbox_inches='tight', dpi=175); plt.close('all')

# Save the cellType information in external file
cellTypesDF = pd.DataFrame(adata.obs['cellType'])
cellTypesDF.to_csv("{0}/03_{1}_cellTypes.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="cellId")


#########################################
# Save session
import dill
filename = "{0}/{1}.pkl".format(output_dir, projName)
dill.dump_session(filename)


# and to load the session again:
import dill
filename = "{0}/{1}.pkl".format(output_dir, projName)
dill.load_session(filename)


















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
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_0_1_2.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.rank_genes_groups(adata, key='rank_genes_r0.5', groups=['3','4','5'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_3_4_5.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.rank_genes_groups(adata, key='rank_genes_r0.5', groups=['6', '7', '8'], fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_cluster_7_6_8.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Pairwise between 0,2,5,10
for l in list(itertools.combinations(['0','2','5','10'], 2)):
  grp = l[0]
  ref = l[1]
  print(list(l))
  sc.tl.rank_genes_groups(adata, groupby='louvain_r1', key_added='rank_genes_{0}_{1}'.format(grp, ref), groups = list(grp), reference = '{}'.format(ref))
  sc.pl.rank_genes_groups(adata, key='rank_genes_{0}_{1}'.format(grp, ref), fontsize=12, show=False)
  plt.savefig("{0}/{1}_marker_genes_ranking_pairwise_{2}_{3}.png".format(plotsDir, bname, grp, ref) , bbox_inches='tight', dpi=175); plt.close('all')




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
plt.savefig("{0}/{1}_marker_genes_cell_annotation_norm_heatmap.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

ma_cell_annotation_norm = sc.tl.marker_gene_overlap(adata, ma_marker_genes, key='rank_genes_r0.5', normalize='reference')
sns.heatmap(ma_cell_annotation_norm, cbar=False, annot=True)
plt.savefig("{0}/{1}_stomach_marker_list_mouse_cellatlas_V1_annotation_norm_heatmap.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

sc.pl.umap(adata, color=['louvain_r1.5', 'Vim', 'Mdm2', 'Trp53', 'Irf8', 'Myc', 'Gamt', 'E2f1', 'Pcna', 'Tgfbr2'], use_raw=False, color_map=mymap)
plt.savefig("{0}/{1}_marker_genes_adult_stomach_mouse_cell_atlas_UMAPs.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Acinar
sc.pl.umap(adata, color=['Ctrb1','Cpa1','Cpb1' ,'Cela2a', 'Cela1', 'Try4','Try5', 'Pnlip', ], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)
sc.pl.umap(adata, color=['Birc5', 'Casp3', 'Stat3', 'Alb'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

sc.pl.umap(adata, color=['Birc5', 'Stat3', 'Reg1', 'Gm26917', 'Ctrb1', 'Clps', 'Hbb-bs', 'Hba-a1','Hba-a2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

sc.pl.umap(adata, color=['Birc5', 'Stat3', 'Reg1', 'Gm26917', 'Ctrb1', 'Clps', 'Hbb-bs', 'Hba-a1','Hba-a2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)
sc.pl.umap(adata, color=['Birc5', 'Blvrb', 'Car2', 'Hbb-bt', 'Clps', 'Hbb-bs', 'Hba-a1','Hba-a2', 'Hspa1b', 'Apoe', 'C1qb', 'Cxcl2', 'Slpi'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)
sc.pl.umap(adata, color=['Sh2d1a', 'Cd3d','Cd3e','Cd8a','Retnlg','S100a8','S100a9','Cxcl2', 'Slpi', 'Srgn', 'Cd84', 'Stip1','Cd44', 'Jak1'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

sc.pl.umap(adata, color=['Btg1', 'Morf4l2','Marcks','Ptprc','BC005537','Ctnnb1','Ptma','AW112010', 'Hnrnpf', 'Hspa1b', 'Hnrnph1', 'Dazap2','Laptm5', 'Id2'], use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9)

sc.pl.umap(adata, color=['louvain_r1.5','Apoe', 'Dcn', 'Cd74','Pf4', 'Lyz2', 'Ly6c1', 'Plvap', 'Fabp5', 'Birc5', 'Ube2c', 'Dmbt1', 'Igfbp3', 'Clps'], use_raw=False, color_map='hot', show=False)
plt.savefig("{0}/{1}_marker_genes_adult_stomach_mouse_cell_atlas_UMAPs.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')


# # Check expression of enterocyte markers
# #Collate all enterocyte markers and get the gene IDs in the data set
# ids_pancreas = np.in1d(adata.var_names, marker_genes['Pancreas'])

# #Calculate the mean expression of enterocyte markers
# adata.obs['Pancreas_marker_expr'] = adata.X[:,ids_pancreas].mean(1)

# #Plot enterocyte expression
# sc.pl.violin(adata, 'Pancreas_marker_expr', groupby='louvain_r0.5')
# sc.pl.umap(adata, color='Pancreas_marker_expr', color_map=mymap)


# On the subclusters
sbadata =adata.copy()
# sc.pl.umap(sbadata, color='louvain', palette=sc.pl.palettes.vega_20, size=50)
plt.figure(figsize=(100,8))
sc.pl.umap(sbadata, color='louvain_r1', palette=sc.pl.palettes.vega_20, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/{1}_Louvainr05_UMAPs.png".format(plotsDir, bname, k) , bbox_inches='tight', dpi=175, dpi=300); plt.close('all')

# Groups interested = 0, 2, 5, 10
sc.tl.rank_genes_groups(sbadata, groupby='louvain_r1', key_added='rank_genes_0_2', groups = '2', reference = '0')
sc.pl.rank_genes_groups(sbadata, key='rank_genes_0_2', fontsize=12)


grp = ['4', '6', '13']
ref = '11'
sc.tl.rank_genes_groups(adata, groupby='louvain_r1', key_added='rank_genes_{0}_{1}'.format(grp, ref), groups = list(grp), reference =  '{}'.format(ref))
sc.pl.rank_genes_groups(adata, key='rank_genes_{0}_{1}'.format(grp, ref), fontsize=12, show=False)
plt.savefig("{0}/{1}_marker_genes_ranking_pairwise_{2}_{3}.png".format(plotsDir, bname, grp, ref) , bbox_inches='tight', dpi=175); plt.close('all')

# Plot umap for the top genes for each cluster
clusterDir = "{0}/clusterDir".format(plotsDir); create_dir(clusterDir)
for cid in [0, 2, 5, 10]:
  cluster = adata[(adata.obs['louvain_r1'] == '{0}'.format(cid))].to_df() 
  topGenes = np.mean(cluster, axis=0).sort_values(ascending=False)[0:25].index.tolist()
  print("\n- Cluster{0}: {1}".format(cid, topGenes))
  plt.figure(figsize=(100,8))
  sc.pl.umap(adata, color=topGenes, use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  plt.savefig("{0}/{1}_top_25_genes_cluster{2}_UMAPs.png".format(clusterDir, bname, cid) , bbox_inches='tight', dpi=175); plt.close('all')
  
# Plot umap for the top marker genes for each cluster
markerDir = "{0}/markerDir/pdf".format(plotsDir); create_dir(markerDir)
for cid in range(len(adata.uns['rank_genes_r1']['names'])):
  topGenes = adata.uns['rank_genes_r1']['names'][cid].tolist()
  print("\n- Group{0}: {1}".format(cid, topGenes))
  # plt.figure(figsize=(100,8))
  # sc.pl.umap(adata, color=topGenes, use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  # plt.savefig("{0}/{1}_marker_genes_group{2}_UMAPs.png".format(markerDir, bname, cid) , bbox_inches='tight', dpi=175); plt.close('all')
  sc.pl.umap(adata, color=topGenes, use_raw=False, color_map=mymap, size=50, legend_loc='on data', edgecolor=None, show=False)
  plt.savefig("{0}/{1}_marker_genes_group{2}_UMAPs.pdf".format(markerDir, bname, cid) , bbox_inches='tight', dpi=175); plt.close('all')
  
# Assign the celltypes to clusters
adata.rename_categories('louvain_r1', ['Birc5⁺/basal', 'Acinar', 'Gastric mucosa/Basal', 'Tcells', 'Dendritic', '5', 'Macrophages', 'Endothelial/Epithelial_Igfbp3⁺', 'Erythrocytes', 'Fibroblasts', '10', 'Restin', ' ', 'Dendritic/Macrophages'])
adata.obs['louvain_r1'].value_counts()
plt.figure(figsize=(100,8))
sc.pl.umap(adata, color='louvain_r1', palette=sc.pl.palettes.vega_20, size=50, legend_loc='on data', edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/{1}_Clusters_CellTypes_UMAPs.png".format(plotsDir, bname, k) , bbox_inches='tight', dpi=175, dpi=300); plt.close('all')


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
filename = "{0}/{1}.pkl".format(output_dir, projName)
dill.dump_session(filename)

# and to load the session again:
import dill
filename = "{0}/{1}.pkl".format(output_dir, projName)
dill.load_session(filename)
#########################################


import scvelo as scv
scv.logging.print_version()

scv.utils.show_proportions(adata)
adata

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)