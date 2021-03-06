# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

# Source: https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb

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
projName        = "bulk1001"
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/individual_samples/{0}".format(projName); create_dir("{0}".format(output_dir))
ccGenes_macosko = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/macosko_cell_cycle_genes_mmu.txt"
ccGenes_regev   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/regev_lab_cell_cycle_genes_mmu.txt"
minGenesPerCell = 5
minCountPerCell = 50
maxCountPerCell = 100000 
minCellsPergene = 2
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
# Running Scanpy 1.4.5.1, on 2020-05-08 02:18.
# scanpy==1.4.5.1 anndata==0.7.1 umap==0.3.10 numpy==1.17.3 scipy==1.4.1 pandas==1.0.3 scikit-learn==0.21.3 statsmodels==0.10.1 python-igraph==0.7.1 louvain==0.6.1

# 1) Reading and performing QC on individual datasets
# 1.1) Reading the data in the anndata object individually
adata = sc.read_10x_h5('input/manec/tissue/bulk1001_mouse_filtered_feature_bc_matrix.h5')

# 1.2) Make the variable names unique and calculate some general qc-stats for genes and cells
adata.var_names_make_unique()

# Convert the sparse count matrices to dense represntation
adata.X = adata.X.toarray()

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

# Total number of cells: 1852
# Number of cells after min count filter: 1852
# Number of cells after max count filter: 1852
# Number of cells after MT filter  : 1845
# Number of cells after Ribo filter: 971
# Number of cells after gene filter: 971
# Total number of genes: 55488
# Number of genes after minCellsPergene filter: 16855
# Number of highly variable genes: 3565

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
ax = fig.add_subplot(2, 3, 6, projection='3d'); sc.pl.umap(adata  , legend_loc=None, ax=ax, color="rb_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="rb_frac UMAP")
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
rawadata = adata.copy()
# 2) # Normalization using SCRAN
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
# Calculations for the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata, random_state = 2105, n_components=3)

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
# normadata = adata.copy()

###############################################################################################################
adata = normadata.copy()

# 5) Plot UMAP for all the tumor cells 
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
plt.savefig("{0}/03_{1}_Tumor_hgMycIresCd2_CellIDs_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# --------------------------------------------------------
# Get tumors for all individual vectors ('humanMyc', 'gap', 'ires', 'humanCd2')
# Unique cell barcode list of list
ucblofl = list()
ucbd    = defaultdict(list)
for t in ['humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2']:
  # Cell barcode data frame
  cbDF = pd.read_csv('/media/rad/HDD2/temp_manec/{1}_hgMycIresCd2_{0}_cellIDs.txt'.format(t,projName), sep="\t", header=None).values.tolist()
  # Unique cell barcode list
  ucbl = get_unique_list(sum(cbDF, []))
  ucbd[t] = ucbl
  ucblofl.append(ucbl)
  print("{0}: {1}".format(t,len(ucbl)))

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

# Set new parameters
size=6
params = {'legend.fontsize': 'large',
          'axes.labelsize' : size,
          'axes.titlesize' : size,
          'xtick.labelsize': size,
          'ytick.labelsize': size,
          'axes.titlepad'  : 15,
          'font.size': size}
plt.rcParams.update(params)
fig = plt.figure(figsize=(10,3))
upsetplot.plot(upsetContent, show_counts='%d', show_percentages=True, fig=fig)
plt.savefig("{0}/03_{1}_Tumor_hgMycIresCd2_individual_CellIDs_upsetPlot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
# Reset to default
mpl.rcParams.update(mpl.rcParamsDefault)

# Plot all hgMycIresCd2 vector components in separate UMAP
fig = plt.figure(figsize=(16,6))
fig.suptitle('hgMycIresCd2 vector components')
# 2D projection
sc.pl.umap(adata, color=['hgMycIresCd2','humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2'], color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, show=False)
plt.savefig("{0}/03_{1}_Tumor_hgMycIresCd2_Components_UMAP_2D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
# 3D projection
sc.pl.umap(adata, color=['hgMycIresCd2','humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2'], color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/03_{1}_Tumor_hgMycIresCd2_Components_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Get tumor cells (and in percentage)
tumorDF = adata.obs[['humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2','hgMycIresCd2']].copy()
tumorDF.to_csv("{0}/05_{1}_tumorComponents_per_vector.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="tissueID")

# 5.5) Save the tumor assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/03_tumorIDs_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/03_tumorIDs_{1}_adata.h5ad" .format(dataDir, projName); tumoridadata  = sc.read_h5ad(adatafile)

#########################################################################
rawadatafile = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName)
rawadata     = sc.read_h5ad(rawadatafile)
adata        = tumoridadata.copy()

# 7) Clustering
# 7.1) Perform clustering - using highly variable genes
sc.tl.louvain(adata, key_added='louvain', random_state=2105)
sc.tl.louvain(adata, resolution=1, key_added='louvain_r1', random_state=2105)
sc.tl.louvain(adata, resolution=1.5, key_added='louvain_r1.5', random_state=2105)
sc.tl.louvain(adata, resolution=2.0, key_added='louvain_r2', random_state=2105)

for i in np.arange(0.1,1,0.1):
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
plt.savefig("{0}/02_norm_{1}_clustering_all_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/02_norm_{1}_clustering_all_louvain_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

cluster_key   = "louvain_r1"
cluster_bname = "louvain_r1"
fig = plt.figure(figsize=(32,12))
# 2D projection
ax = fig.add_subplot(2, 4, 1);                  sc.pl.umap(rawadata,                  ax=ax, color="log_counts", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw log_counts UMAP")
ax = fig.add_subplot(2, 4, 2);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="log_counts", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Norm log_counts UMAP")
ax = fig.add_subplot(2, 4, 3);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Norm mt_frac UMAP")
ax = fig.add_subplot(2, 4, 4);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="rb_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Norm rb_frac UMAP")
# 3D projection
ax = fig.add_subplot(2, 4, 5, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="log_counts", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw log_counts UMAP")
ax = fig.add_subplot(2, 4, 6, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color="log_counts", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Norm log_counts UMAP")
ax = fig.add_subplot(2, 4, 7, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Norm mt_frac UMAP")
ax = fig.add_subplot(2, 4, 8, projection='3d'); sc.pl.umap(adata  , legend_loc=None, ax=ax, color="rb_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Norm rb_frac UMAP")
plt.tight_layout()
plt.savefig("{0}/02_norm_{1}_{2}_tissueID_counts_mtfrac_UMAP.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Louvain UMAPs
fig = plt.figure(figsize=(16,6))
fig.suptitle("{0} UMAP".format(cluster_key))
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  
sc.pl.umap(adata, legend_loc=None, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.5, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(1, 2, 2, projection='3d'); 
sc.pl.umap(adata, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.5, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/02_norm_{1}_clustering_{2}_UMAP_2D3D.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
cluster_key_groups = adata.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = adata.obs[cluster_key].value_counts().to_dict()

# Louvain UMAPs
subplot_title_fontsize = 12
subplot_title_width    = 50
ncols  = len(cluster_key_groups) + 1
fig = plt.figure(figsize=(14, 7*ncols))
fig.suptitle("{0} UMAP".format(cluster_key))
# Main Louvain Cluster
ax = fig.add_subplot(ncols,2, 1); sc.pl.umap(adata, legend_loc=None, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(ncols,2, 2, projection='3d'); sc.pl.umap(adata, ax=ax, color=cluster_key, palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
# Partial visualizaton of a subset of groups in embedding
m=3; n=4
for i,b in enumerate(cluster_key_groups):
  print(i, b)
  ax = fig.add_subplot(ncols,2, i+m);                  sc.pl.umap(adata, legend_loc=None, ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
  ax = fig.add_subplot(ncols,2, i+n, projection='3d'); sc.pl.umap(adata                 , ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
  m+=1; n+=1

plt.tight_layout()
plt.savefig("{0}/02_norm_{1}_clustering_{2}_UMAP_individual_clusters.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

##################################################################
# 8.1) Marker genes & cluster annotation
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby=cluster_key, key_added='rank_genes_{0}'.format(cluster_key), n_genes=adata.shape[1])

# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_{0}'.format(cluster_key), fontsize=12, show=False)
plt.savefig("{0}/02_norm_{1}_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# Annotation of cluster r_0.5 with known marker genes
markerDir = "{0}/markerDir".format(plotsDir); create_dir(markerDir)
subplot_title_fontsize = 12
subplot_title_width    = 50

# Read the marker genes into a pandas dataframe
# marker_file  = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_V2.txt'
# markersDF    = pd.read_csv(marker_file, sep="\t")
# marker_genes = markersDF.groupby('CellLines')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()
# marker_genes_cellTypes = markersDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()

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
  plt.savefig("{0}/31_{1}_marker_genes_stomach_{2}_UMAPs.png".format(markerDir, bname, k) , bbox_inches='tight', dpi=100); plt.close('all')

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

# Get clean marker dict
adata_expressed_genes = adata.var.index.tolist()
marker_genes_filtered_dict = defaultdict()
for k,v in marker_genes_cellTypes.items():
  new_genes_list = [x for x in v if x in adata_expressed_genes]
  if new_genes_list:
    marker_genes_filtered_dict[k] = new_genes_list

ma_marker_genes_filtered_dict = defaultdict()
for k,v in ma_marker_genes.items():
  new_genes_list = [x for x in v if x in adata_expressed_genes]
  if new_genes_list:
    ma_marker_genes_filtered_dict[k] = new_genes_list

# 8.2) Other marker gene visualization
marker_list_name = "stomach_V2"
# 8.2.1) Dot plots: The dotplot visualization provides a compact way of showing per group, the fraction of cells expressing a gene (dot size) and the mean expression of the gene in those cell (color scale).
# The use of the dotplot is only meaningful when the counts matrix contains zeros representing no gene counts. dotplot visualization does not work for scaled or corrected matrices in which cero counts had been replaced by other values.
sc.pl.dotplot(adata, marker_genes_filtered_dict, groupby=cluster_key, log=True, figsize=(40,12), show=False, dendrogram=True)
plt.savefig("{0}/02_norm_{1}_{2}_31_marker_genes_{3}_dotplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.2) Matrix plots: The matrixplot shows the mean expression of a gene in a group by category as a heatmap. In contrast to dotplot, the matrix plot can be used with corrected and/or scaled counts. By default raw counts are used.
sc.pl.matrixplot(adata, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False,cmap='Reds',  figsize=(40,12), standard_scale='group', show=False)
plt.savefig("{0}/02_norm_{1}_{2}_31_marker_genes_{3}_scaled_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.matrixplot(adata, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False, cmap='Reds', figsize=(40,12), standard_scale='group', vmin=0.5, show=False)
plt.savefig("{0}/02_norm_{1}_{2}_31_marker_genes_{3}_scaled_vmin0_05_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.3) Tracksplots: The track plot shows the same information as the heatmap, but, instead of a color scale, the gene expression is represented by height.
ad = adata.copy()
ad.raw.X.data = np.exp(ad.raw.X.data)
ax = sc.pl.tracksplot(ad, marker_genes_filtered_dict, groupby=cluster_key, log=True, dendrogram=True, show=False, figsize=(50,30))
plt.savefig("{0}/02_norm_{1}_{2}_31_marker_genes_{3}_tracksplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')

marker_list_name = "mouse_cellatlas"
# 8.2.1) Dot plots
sc.pl.dotplot(adata, ma_marker_genes_filtered_dict, groupby=cluster_key, log=True, figsize=(40,12), show=False, dendrogram=True)
plt.savefig("{0}/02_norm_{1}_{2}_32_marker_genes_{3}_dotplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.2) Matrix plots
sc.pl.matrixplot(adata, ma_marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False,cmap='Reds',  figsize=(40,12), standard_scale='group', show=False)
plt.savefig("{0}/02_norm_{1}_{2}_32_marker_genes_{3}_scaled_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.matrixplot(adata, ma_marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False, cmap='Reds', figsize=(40,12), standard_scale='group', vmin=0.5, show=False)
plt.savefig("{0}/02_norm_{1}_{2}_32_marker_genes_{3}_scaled_vmin0_05_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.3) Tracksplots
ax = sc.pl.tracksplot(ad, ma_marker_genes_filtered_dict, groupby=cluster_key, log=True, dendrogram=True, show=False, figsize=(50,30))
plt.savefig("{0}/02_norm_{1}_{2}_32_marker_genes_{3}_tracksplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')

# 8.4) Dataframe of ranked genes
# Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
cluster_key        = "louvain_r1"
cluster_bname      = "louvain_r1"
cluster_key_groups = adata.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = adata.obs[cluster_key].value_counts().to_dict()
rankGenesDir       = "{0}/rankedGenes/{1}".format(dataDir,cluster_bname); create_dir(rankGenesDir)
for g in cluster_key_groups:
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(adata.uns['rank_genes_{0}'.format(cluster_key)][n])[g]
  # Save dataframes
  ngDF.to_csv("{0}/03_{1}_rank_genes_{2}_{3}.txt".format(rankGenesDir, projName, cluster_bname, g), sep='\t', header=True, index=False, float_format='%.2g')

# 8.6) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); markeradata  = sc.read_h5ad(adatafile)
# adata = markeradata.copy()
# Finished on 2020-04Apr-27 00:33
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
