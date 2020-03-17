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
projName        = "manec_tissues_merged_except1079" # MANEC_merged_except1079_hMYC_forcecells
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/01_tissue/{0}".format(projName); create_dir("{0}".format(output_dir))
cc_genes_file   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/Macosko_cell_cycle_genes.txt"
minGenesPerCell = 100
minCountPerCell = 300
minCellsPergene = 25
bname           = projName
topqcDir        = "{0}/analysis".format(output_dir); create_dir(topqcDir)
topcountsDir    = "{0}/data".format(output_dir)    ; create_dir(topcountsDir)

# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

# 1 Reading the data
# https://github.com/theislab/scanpy/issues/267
tissueFilenames = [
                    'input/manec/tissue/bulk997_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/bulk1001_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/bulk1018_mouse_filtered_feature_bc_matrix.h5', 
                    'input/manec/tissue/stomach1001_mouse_filtered_feature_bc_matrix.h5'
                  ]
adatas          = [sc.read_10x_h5(f) for f in tissueFilenames]

# Get the dictionary of tissue ids
tissueIdDict =  {
                    '0':'bulk1001', 
                    '1':'bulk997', 
                    '2':'bulk1018', 
                    '3':'stomach1001'
                }

# Make the variable names unique for each annData object separately and run DCA on them
for i, adata in enumerate(adatas):
  qcDir           = "{0}/analysis".format(output_dir); create_dir(qcDir)

  # I will run in this in the loop at the end
  i         = 2
  adata     = adatas[2]
  qcDir     = "{0}/{1}".format(topqcDir, tid)    ; create_dir(qcDir)
  countsDir = "{0}/{1}".format(topcountsDir, tid); create_dir(countsDir)

  tid = tissueIdDict[str(i)]
  adata.var_names_make_unique()

  # Perform basic prefiltering
  print("- {0} with shape {1}".format(tid, adata.to_df().shape))

  # Convert the sparse count matrices to dense represntation
  adata.X = adata.X.toarray()

  # AnnData object with n_obs × n_vars = 9038 × 31053 
  #     obs: 'batch', 'tissueID'
  #     var: 'gene_ids', 'feature_types', 'genome'

  # Checking the total size of the data set
  adata.shape # We have 9038 cells and 31053 genes in the dataset

  # 2) Quality control 
  # 2.1) Calculate QC covariates
  adata.obs['n_counts'] = adata.X.sum(1)
  adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
  adata.obs['n_genes'] = (adata.X > 0).sum(1)

  mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
  adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts']

  # 2.2) Data quality summary plots
  p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', show=False)
  plt.savefig("{0}/01_raw_{1}_genes_counts_scatterplot.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

  p2 = sc.pl.scatter(adata[adata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', show=False)
  plt.savefig("{0}/01_raw_{1}_genes_counts_scatterplot_zoomedin.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

  # 2.3) Thresholding decision based on counts
  p3 = sns.distplot(adata.obs['n_counts'], kde=False); #plt.show()
  plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  p4 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<1000], kde=False, bins=100); #plt.show()
  plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot_lessthan_1000.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  p5 = sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>1000], kde=False, bins=100); #plt.show()
  plt.savefig("{0}/01_raw_{1}_ncounts_histogramplot_greaterthan_1000.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

  # 2.4) Thresholding decision based on genes
  p6 = sns.distplot(adata.obs['n_genes'], kde=False, bins=100); # plt.show()
  plt.savefig("{0}/01_raw_{1}_genes_histogramplot.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  p7 = sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<500], kde=False, bins=100); # plt.show()
  plt.savefig("{0}/01_raw_{1}_genes_histogramplot_lessthan_500.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

  # 2.5) Filter cells according to identified QC thresholds:
  origadata = adata.copy()
  print('Total number of cells: {:d}'.format(adata.n_obs))

  # Reset the thresholds
  minGenesPerCell = 350
  minCountPerCell = 300
  minCellsPergene = 25

  sc.pp.filter_cells(adata, min_counts = minCountPerCell)
  print('Number of cells after min count filter: {:d}'.format(adata.n_obs))

  # sc.pp.filter_cells(adata, max_counts = 40000)
  # print('Number of cells after max count filter: {:d}'.format(adata.n_obs))

  adata = adata[adata.obs['mt_frac'] < 0.25]
  print('Number of cells after MT filter: {:d}'.format(adata.n_obs))

  sc.pp.filter_cells(adata, min_genes = minGenesPerCell)
  print('Number of cells after gene filter: {:d}'.format(adata.n_obs))

  # bulk1018
  # Total number of cells: 897
  # Number of cells after min count filter: 897
  # Number of cells after MT filter: 818
  # Trying to set attribute `.obs` of view, copying.
  # Number of cells after gene filter: 688

  # 2.6) Filter genes according to identified QC thresholds:
  # Min 5 cells - filters out 0 count genes
  print('Total number of genes: {:d}'.format(adata.n_vars))
  sc.pp.filter_genes(adata, min_cells=minCellsPergene)
  print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
  # Total number of genes: 11664
  # Number of genes after cell filter: 8068

  # 2.8) Calculations for the visualizations
  sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
  print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
  sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
  sc.pp.neighbors(adata, random_state = 2105)
  sc.tl.umap(adata, random_state = 2105, n_components=3)

  # 2.9) Plot visualizations
  sc.pl.pca_scatter(adata, color='n_counts',show=False)
  plt.savefig("{0}/01_raw_{1}_clustering_ncounts_PCA.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  sc.pl.umap(adata, color=['mt_frac'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  plt.savefig("{0}/01_raw_{1}_clustering_mt_frac_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  sc.pl.umap(adata, color=['mt_frac'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
  plt.savefig("{0}/01_raw_{1}_clustering_mt_frac_UMAP_3D.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

  p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac', show=False)
  plt.savefig("{0}/01.2_raw_{1}_genes_counts_scatterplot.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

  # 3) Normalization using SCRAN
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
  plt.savefig("{0}/02_norm_{1}_scrna_sizefactors_vs_ncounts.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  sc.pl.scatter(adata, 'size_factors', 'n_genes' , show=False)
  plt.savefig("{0}/02_norm_{1}_scrna_sizefactors_vs_ngenes.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  sns.distplot(size_factors, bins=50, kde=False)
  plt.savefig("{0}/02_norm_{1}_scrna_sizefactors_histogram.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

  # Keep the count data in a counts layer
  adata.layers["counts"] = adata.X.copy()

  # Normalize adata 
  adata.X /= adata.obs['size_factors'].values[:,None]
  sc.pp.log1p(adata)

  # Store the full data set in 'raw' as log-normalised data for statistical testing
  adata.raw = adata

  # Highly Variable Genes
  sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
  print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
  sc.pl.highly_variable_genes(adata, show=False)
  plt.savefig("{0}/02_norm_{1}_highly_variable_genes.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

  # 4) Biological correction: Cell cycle scoring
  # Score cell cycle and visualize the effect
  cc_genes         = pd.read_table(cc_genes_file, delimiter='\t')
  s_genes          = cc_genes['S'].dropna()
  g2m_genes        = cc_genes['G2.M'].dropna()
  s_genes_mm       = [gene.lower().capitalize() for gene in s_genes]
  g2m_genes_mm     = [gene.lower().capitalize() for gene in g2m_genes]
  s_genes_mm_ens   = adata.var_names[np.in1d(adata.var_names, s_genes_mm)]
  g2m_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, g2m_genes_mm)]

  sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)
  sc.pl.umap(adata, color=['S_score', 'G2M_score'], use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  plt.savefig("{0}/02_norm_{1}_S_G2M_Phase_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  sc.pl.umap(adata, color='phase', use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
  plt.savefig("{0}/02_norm_{1}_combined_Phase_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  sc.pl.umap(adata, color='phase', use_raw=False, palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
  plt.savefig("{0}/02_norm_{1}_combined_Phase_UMAP_3D.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

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

  # 0    563
  # 1    125
  # Name: louvain_r0.1, dtype: int64

  # 0    487
  # 1    125
  # 2     76
  # Name: louvain_r0.2, dtype: int64

  # 0    297
  # 1    184
  # 2    125
  # 3     82
  # Name: louvain_r0.3, dtype: int64

  # 0    280
  # 1    193
  # 2    125
  # 3     90
  # Name: louvain_r0.4, dtype: int64

  # 0    194
  # 1    144
  # 2    134
  # 3    125
  # 4     91
  # Name: louvain_r0.5, dtype: int64

  # 0    141
  # 1    134
  # 2    125
  # 3    114
  # 4     93
  # 5     81
  # Name: louvain_r0.6, dtype: int64

  # 0    168
  # 1    125
  # 2    118
  # 3     93
  # 4     68
  # 5     65
  # 6     51
  # Name: louvain_r0.7, dtype: int64

  # 0    148
  # 1    125
  # 2    120
  # 3     95
  # 4     68
  # 5     67
  # 6     65
  # Name: louvain_r0.8, dtype: int64

  # 0    148
  # 1    113
  # 2     97
  # 3     95
  # 4     76
  # 5     66
  # 6     65
  # 7     28
  # Name: louvain_r0.9, dtype: int64


  print(adata.obs['louvain'].value_counts())
  print(adata.obs['louvain_r1'].value_counts())
  print(adata.obs['louvain_r1.5'].value_counts())
  print(adata.obs['louvain_r2'].value_counts())

  # 0    144
  # 1     97
  # 2     96
  # 3     94
  # 4     76
  # 5     66
  # 6     62
  # 7     28
  # 8     25
  # Name: louvain, dtype: int64

  # 0    144
  # 1     97
  # 2     96
  # 3     94
  # 4     76
  # 5     66
  # 6     62
  # 7     28
  # 8     25
  # Name: louvain_r1, dtype: int64

  # 0    119
  # 1     98
  # 2     95
  # 3     83
  # 4     76
  # 5     68
  # 6     62
  # 7     33
  # 8     28
  # 9     26
  # Name: louvain_r1.5, dtype: int64

  # 0     84
  # 1     81
  # 3     63
  # 2     63
  # 4     62
  # 5     59
  # 6     57
  # 7     49
  # 9     39
  # 8     39
  # 10    37
  # 11    28
  # 12    27
  # Name: louvain_r2, dtype: int64


# Calculations for the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.umap(adata, random_state = 2105, n_components=3)

# Plot visualizations
# Visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=['louvain'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_louvain_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_all_louvain_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Looking at the data, choose r = 0.5
sc.pl.umap(adata, color=['louvain_r0.5'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/02_norm_{1}_clustering_louvain_r05_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

sc.pl.umap(adata, color=['louvain_r0.5'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, legend_loc='on data', show=False)
plt.savefig("{0}/02_norm_{1}_clustering_louvain_r05_legend_onData_UMAP.png".format(qcDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

