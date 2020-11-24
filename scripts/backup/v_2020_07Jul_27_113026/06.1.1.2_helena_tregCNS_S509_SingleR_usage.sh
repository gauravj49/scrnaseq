# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq
# Source: https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/Case-study_Mouse-intestinal-epithelium_1906.ipynb
# NOTE: 
# 1) plt.subplot(number_of_rows, number_of_columns, index_of_the_subplot) 

##################################################################
ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)

# Loading the python libraries environment
%load scripts/load_python_modules.py
%load scripts/scrnaseq_module.py

# System variables and directories
projName        = "S509"
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/tregCNS/individual_samples/{0}".format(projName); create_dir("{0}".format(output_dir))
ccGenes_macosko = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/macosko_cell_cycle_genes_mmu.txt"
ccGenes_regev   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/regev_lab_cell_cycle_genes_mmu.txt"
minGenesPerCell = 5
minCountPerCell = 50
maxCountPerCell = 100000 
minCellsPergene = 2
mtGenesFilter   = 0.25
rbGenesFilter   = 0.30
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
# Memory usage: current 1.07 GB, difference +1.07 GB
# Running Scanpy 1.4.6, on 2020-07-13 10:35.
# scanpy==1.4.6 anndata==0.7.3 umap==0.4.4 numpy==1.19.0rc2 scipy==1.4.1 pandas==0.25.3 scikit-learn==0.23.1 statsmodels==0.11.1

#############################################################################
# Working after Louvain clustering
# 8.6) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
# adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# Read back the corrected adata object
adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); markeradata  = sc.read_h5ad(adatafile)
adata = markeradata.copy()

# adata = markeradata.copy()
cluster_key        = "louvain_r1"
cluster_bname      = "louvain_r1"
cluster_key_groups = adata.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = adata.obs[cluster_key].value_counts().to_dict()














































# 8) CellType assignments
# 8.1) Add new categories
# Get a new subcluster column
# # To exclude
#  1 = Macrophages
#  5 = Erythrocytes
#  6 = Endothelial_Epithelial_Igfbp3pos
#  7 = Fibroblasts
##### to keep for subgroup #####
#  0 = Cluster0
#  2 = Cluster2
#  3 = Tumor3
#  4 = Tumor4 
#  8 = ECL
# Categories to rename
adata.obs[cluster_key].cat.categories
# Get a new cell type column from the annotation of the louvain_r0.5 clusters
adata.obs['cellType'] = adata.obs[cluster_key]
# Add new categories
adata.obs['cellType'].cat.add_categories(['Cluster0','Macrophages','Cluster2','Tumor3','Tumor4','Erythrocytes','Endothelial_Epithelial_Igfbp3pos','Fibroblasts','ECL'], inplace=True) 

adata.obs['cellType'].loc[adata.obs['cellType' ]=='0' ]  = 'Cluster0'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='1' ]  = 'Macrophages'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='2' ]  = 'Cluster2'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='3' ]  = 'Tumor3'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='4' ]  = 'Tumor4'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='5' ]  = 'Erythrocytes'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='6' ]  = 'Endothelial_Epithelial_Igfbp3pos'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='7' ]  = 'Fibroblasts'
adata.obs['cellType'].loc[adata.obs['cellType' ]=='8' ]  = 'ECL'

# Remove old categories
adata.obs['cellType'].cat.remove_categories(adata.obs[cluster_key].cat.categories.tolist(), inplace=True)
# List new categories
adata.obs['cellType'].cat.categories

# Draw Umaps with the categories
fig = plt.figure(figsize=(36,20))
fig.suptitle('CellType UMAP')
# 2D projection
ax = fig.add_subplot(2, 3, 1);                  sc.pl.umap(adata, legend_loc='on data', ax=ax, color="cellType", palette=sc.pl.palettes.vega_20, size=150, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(2, 3, 2);                  sc.pl.umap(adata, legend_loc=None     , ax=ax, color="cellType",  palette=sc.pl.palettes.vega_20, size=150, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(2, 3, 3, projection='3d'); sc.pl.umap(adata                      , ax=ax, color="cellType", palette=sc.pl.palettes.vega_20, size=150, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
# Save the UMAP
plt.savefig("{0}/04_{1}_clustering_CellType_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=200); plt.close('all')

# 8.3) Calculate marker genes (one vs. rest)
sc.tl.rank_genes_groups(adata, groupby='cellType', key_added='rank_genes_cellTypes', n_genes=adata.shape[1])
# Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_cellTypes', fontsize=12, show=False)
plt.savefig("{0}/04_subGroup_{1}_{2}_marker_genes_ranking_cellType.png".format(plotsDir, bname, 'cellType') , bbox_inches='tight', dpi=175); plt.close('all')

# Get average expression DF
# In [23]: np.mean(adata.to_df(), axis=0)
# Out[23]:
# Sox17      0.039553
# Mrpl15     0.152230
#              ...
# mt-Nd4l    0.155424
# mt-Nd4     1.708211
# Length: 8473, dtype: float32
avgExpDF = pd.DataFrame(np.mean(adata.to_df(), axis=0))
avgExpDF.rename(columns={0: "MeanExpression"}, inplace=True)
# Get all cellTypes into the list
cellTypeCategories = adata.obs['cellType'].cat.categories.tolist()
# Get data dir
rgDataDir  = "{0}/rankedGenes/SubCategories".format(dataDir); create_dir(rgDataDir)
for grp in cellTypeCategories:
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(adata.uns['rank_genes_cellTypes'][n])[grp]
  # Add treatment and reference group name
  ngDF['Treatment'] = grp
  ngDF['Reference'] = 'rest'
  # Convert genes columns to index
  ngDF.set_index('names', inplace=True)
  ngDF['MeanExpression'] = avgExpDF['MeanExpression']
  # Rearragnge columns
  ngDF = ngDF[['MeanExpression', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj','Treatment','Reference']]
  # Save the dataframe
  ngDF.to_csv("{0}/04_{1}_subGroup_{2}_marker_genes_ranking_cellType.txt".format(rgDataDir, projName, grp), sep='\t', header=True, index=True, float_format='%.2g')

# 8.4) Calculate pairwise marker genes list
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
#  1 = Macrophages
#  5 = Erythrocytes
#  6 = Endothelial_Epithelial_Igfbp3pos
#  7 = Fibroblasts
adataSubGroup = adata[~((adata.obs['cellType']=='Macrophages') | (adata.obs['cellType']=='Erythrocytes') |(adata.obs['cellType']=='Endothelial_Epithelial_Igfbp3pos')|(adata.obs['cellType']=='Fibroblasts'))].copy()
# adataSubGroup.shape # (529, 17263)

# Calculations for the visualizations
sc.pp.neighbors(adataSubGroup, random_state = 2105)
sc.tl.umap(adataSubGroup, random_state = 2105, n_components=3)

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
cellBarCodes = pd.read_csv('/media/rad/HDD2/temp_manec/hgMycIresCd2_cellIDs.txt', sep="\t", header=None).values.tolist()
cl  = sum(cellBarCodes, [])
ucl = get_unique_list(cl)
mylist = adataSubGroup.obs.index.values
humaniresmyc = list()
for e in mylist: 
  flag = 0
  for s in ucl: 
      if s in e: 
          flag = 1 
          break
  humaniresmyc.append(flag)

adataSubGroup.obs['hgMycIresCd2'] = humaniresmyc
fig = plt.figure(figsize=(16,6))
fig.suptitle('hgMycIresCd2')
ax = fig.add_subplot(1, 2, 1); sc.pl.umap(adataSubGroup, legend_loc=None, ax=ax, color="hgMycIresCd2", color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(1, 2, 2, projection='3d'); sc.pl.umap(adataSubGroup, ax=ax, color="hgMycIresCd2", color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/05_subGroup_{1}_Tumor_hgMycIresCd2_CellIDs_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Subcluster keys
cluster_key            = "louvain_r1"
cluster_bname          = "louvain_r1"
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

# 9.4) Marker genes & cluster annotation
# Calculate marker genes
sc.tl.rank_genes_groups(adataSubGroup, groupby=cluster_key, key_added='rank_genes_{0}'.format(cluster_key), n_genes=adataSubGroup.shape[1])

# Plot marker genes
sc.pl.rank_genes_groups(adataSubGroup, key='rank_genes_{0}'.format(cluster_key), fontsize=12, show=False)
plt.savefig("{0}/05_subGroup_{1}_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# Annotation of cluster r_0.5 with known marker genes
markerDir = "{0}/markerDir/SubCategories".format(plotsDir); create_dir(markerDir)
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

# Get clean marker dict
adata_expressed_genes = adataSubGroup.var.index.tolist()
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
sc.pl.dotplot(adataSubGroup, marker_genes_filtered_dict, groupby=cluster_key, log=True, figsize=(40,12), show=False, dendrogram=True)
plt.savefig("{0}/05_subGroup_{1}_{2}_31_marker_genes_{3}_dotplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.2) Matrix plots: The matrixplot shows the mean expression of a gene in a group by category as a heatmap. In contrast to dotplot, the matrix plot can be used with corrected and/or scaled counts. By default raw counts are used.
sc.pl.matrixplot(adataSubGroup, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False,cmap='Reds',  figsize=(40,12), standard_scale='group', show=False)
plt.savefig("{0}/05_subGroup_{1}_{2}_31_marker_genes_{3}_scaled_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.matrixplot(adataSubGroup, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False, cmap='Reds', figsize=(40,12), standard_scale='group', vmin=0.5, show=False)
plt.savefig("{0}/05_subGroup_{1}_{2}_31_marker_genes_{3}_scaled_vmin0_05_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.3) Tracksplots: The track plot shows the same information as the heatmap, but, instead of a color scale, the gene expression is represented by height.
ad = adataSubGroup.copy()
ad.raw.X.data = np.exp(ad.raw.X.data)
ax = sc.pl.tracksplot(ad, marker_genes_filtered_dict, groupby=cluster_key, log=True, dendrogram=True, show=False, figsize=(50,30))
plt.savefig("{0}/05_subGroup_{1}_{2}_31_marker_genes_{3}_tracksplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')

marker_list_name = "mouse_cellatlas"
# 8.2.1) Dot plots
sc.pl.dotplot(adataSubGroup, ma_marker_genes_filtered_dict, groupby=cluster_key, log=True, figsize=(40,12), show=False, dendrogram=True)
plt.savefig("{0}/05_subGroup_{1}_{2}_32_marker_genes_{3}_dotplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.2) Matrix plots
sc.pl.matrixplot(adataSubGroup, ma_marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False,cmap='Reds',  figsize=(40,12), standard_scale='group', show=False)
plt.savefig("{0}/05_subGroup_{1}_{2}_32_marker_genes_{3}_scaled_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.matrixplot(adataSubGroup, ma_marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False, cmap='Reds', figsize=(40,12), standard_scale='group', vmin=0.5, show=False)
plt.savefig("{0}/05_subGroup_{1}_{2}_32_marker_genes_{3}_scaled_vmin0_05_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.3) Tracksplots
ax = sc.pl.tracksplot(ad, ma_marker_genes_filtered_dict, groupby=cluster_key, log=True, dendrogram=True, show=False, figsize=(50,30))
plt.savefig("{0}/05_subGroup_{1}_{2}_32_marker_genes_{3}_tracksplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')

# 8.4) Dataframe of ranked genes
# Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
cluster_key        = "louvain_r1"
cluster_bname      = "louvain_r1"
cluster_key_groups = adataSubGroup.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = adataSubGroup.obs[cluster_key].value_counts().to_dict()
rankGenesDir       = "{0}/rankedGenes/SubCategories/{1}".format(dataDir,cluster_bname); create_dir(rankGenesDir)
for g in cluster_key_groups:
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(adataSubGroup.uns['rank_genes_{0}'.format(cluster_key)][n])[g]
  # Save dataframes
  ngDF.to_csv("{0}/04_subGroup_{1}_rank_genes_{2}_{3}.txt".format(rankGenesDir, projName, cluster_bname, g), sep='\t', header=True, index=False, float_format='%.2g')

# 9.5) Save the cellType assigned adataSubGroup into a file
# Write the adataSubGroup and cadataSubGroup object to file
adataSubGroupfile  = "{0}/06_markerGenes_{1}_adataSubGroup.h5ad" .format(dataDir, projName); adataSubGroup.write(adataSubGroupfile)
# # Read back the corrected adataSubGroup object
# adataSubGroupfile  = "{0}/06_markerGenes_{1}_adataSubGroup.h5ad" .format(dataDir, projName); markeradataSubGroup  = sc.read_h5ad(adataSubGroupfile)

################################################################################
# adataSubGroup = markeradataSubGroup.copy()

# 10) subCellType assignments
# Subcluster keys
cluster_key        = "louvain_r1"
cluster_bname      = "louvain_r1"
cluster_key_groups = adataSubGroup.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = adataSubGroup.obs[cluster_key].value_counts().to_dict()
subplot_title_fontsize = 12
subplot_title_width    = 50

# 10.1) Add new categories 
# Categories to rename
adataSubGroup.obs[cluster_key].cat.categories
# Get a new cell type column from the annotation of the louvain_r0.5 clusters
adataSubGroup.obs['subCellType'] = adataSubGroup.obs[cluster_key]
# Add new categories
adataSubGroup.obs['subCellType'].cat.add_categories(['C0_SubCluster0','C1_Macrophages','C2_Tumor2_PitCell','C3_Tumor3','C4_SubCluster4','C5_ECL','C6_SubCluster6'], inplace=True) 
# Get a new subcluster column
# 0 = 'SubCluster0'
# 1 = 'Macrophages'
# 2 = 'Tumor2_PitCell'
# 3 = 'Tumor3'
# 4 = 'SubCluster4'
# 5 = 'ECL'
# 6 = 'SubCluster6'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='0']  = 'C0_SubCluster0'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='1']  = 'C1_Macrophages'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='2']  = 'C2_Tumor2_PitCell'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='3']  = 'C3_Tumor3'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='4']  = 'C4_SubCluster4'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='5']  = 'C5_ECL'
adataSubGroup.obs['subCellType'].loc[adataSubGroup.obs['subCellType' ]=='6']  = 'C6_SubCluster6'

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

# 10.3) Calculate marker genes (one vs. rest)
sc.tl.rank_genes_groups(adataSubGroup, groupby='subCellType', key_added='rank_genes_{0}'.format(cluster_key), n_genes=adataSubGroup.shape[1])
# Plot marker genes
sc.pl.rank_genes_groups(adataSubGroup, key='rank_genes_{0}'.format(cluster_key), fontsize=12, show=False)
plt.savefig("{0}/06_{1}_{2}_marker_genes_ranking_subCellType.png".format(plotsDir, bname, 'subCellType') , bbox_inches='tight', dpi=175); plt.close('all')

# Get average expression DF
# In [23]: np.mean(adata.to_df(), axis=0)
# Out[23]:
# Sox17      0.039553
# Mrpl15     0.152230
#              ...
# mt-Nd4l    0.155424
# mt-Nd4     1.708211
# Length: 8473, dtype: float32
avgExpDF = pd.DataFrame(np.mean(adataSubGroup.to_df(), axis=0))
avgExpDF.rename(columns={0: "MeanExpression"}, inplace=True)
# Get all cellTypes into the list
cellTypeCategories = adataSubGroup.obs['subCellType'].cat.categories.tolist()
# Get data dir
rgDataDir  = "{0}/rankedGenes/SubCategories_oneVsRest".format(dataDir); create_dir(rgDataDir)
for grp in cellTypeCategories:
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(adataSubGroup.uns['rank_genes_{0}'.format(cluster_key)][n])[grp]
  # Add treatment and reference group name
  ngDF['Treatment'] = grp
  ngDF['Reference'] = 'rest'
  # Convert genes columns to index
  ngDF.set_index('names', inplace=True)
  ngDF['MeanExpression'] = avgExpDF['MeanExpression']
  # Rearragnge columns
  ngDF = ngDF[['MeanExpression', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj','Treatment','Reference']]
  # Save the dataframe
  ngDF.to_csv("{0}/04_{1}_subGroup_{2}_marker_genes_ranking_cellType.txt".format(rgDataDir, projName, grp), sep='\t', header=True, index=True, float_format='%.2g')

# 10.4) Calculate pairwise marker genes list
# Get all subCellTypes into the list
subCellTypeCategories = adataSubGroup.obs['subCellType'].cat.categories.tolist()
# Get the list of all unique pairwise combinations
subCellTypePairComb = [comb for comb in combinations(subCellTypeCategories, 2)]
# Get pairwise plots dir
pwsgPlotsDir = "{0}/rankedGenes/SubCategories_pairwise_rankedGenes".format(plotsDir); create_dir(pwsgPlotsDir)
# Get pairwise data dir
pwsgDataDir  = "{0}/SubCategories_pairwise_rankedGenes".format(dataDir); create_dir(pwsgDataDir)
pwlogfcDF    = pd.DataFrame()
pwpvadjDF    = pd.DataFrame()
pwnamesDF    = pd.DataFrame()
# Calculate pairwise marker genes  
for grp,ref in subCellTypePairComb:
  print("- Calculating pairwise marker genes for group_v_reference: {0}_v_{1}".format(grp, ref))
  # Get genes ranking
  keyName = 'rank_genes_{0}_v_{1}'.format(grp,ref)
  sc.tl.rank_genes_groups(adataSubGroup, groupby='subCellType', groups=[grp], key_added=keyName, reference=ref, n_genes=adataSubGroup.shape[1])
  # Get the individual dataframes and add it as a column
  pwnamesDF["{0}_v_{1}".format(grp, ref)] = pd.DataFrame(adataSubGroup.uns[keyName]['names']).values.flatten()
  pwlogfcDF["{0}_v_{1}".format(grp, ref)] = pd.DataFrame(adataSubGroup.uns[keyName]['logfoldchanges']).values.flatten()
  pwpvadjDF["{0}_v_{1}".format(grp, ref)] = pd.DataFrame(adataSubGroup.uns[keyName]['pvals_adj']).values.flatten()
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
  # Convert genes columns to index
  ngDF.set_index('names', inplace=True)
  ngDF['MeanExpression'] = avgExpDF['MeanExpression']
  # Rearragnge columns
  ngDF = ngDF[['MeanExpression', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj','Treatment','Reference']]
  # Save the dataframe
  ngDF.to_csv("{0}/04_{1}_subGroup_{2}_marker_genes_ranking_cellType.txt".format(pwsgDataDir, projName, keyName), sep='\t', header=True, index=True, float_format='%.2g')
  # Add column to logfoldchange DF
# Save the dataframes
pwnamesDF.to_csv("{0}/04_{1}_subGroup_marker_genes_ranking_names.txt".format(pwsgDataDir, projName), sep='\t', header=True, index=True, float_format='%.2g')
pwlogfcDF.to_csv("{0}/04_{1}_subGroup_marker_genes_ranking_logfoldchanges.txt".format(pwsgDataDir, projName), sep='\t', header=True, index=True, float_format='%.2g')
pwpvadjDF.to_csv("{0}/04_{1}_subGroup_marker_genes_ranking_pvals_adj.txt".format(pwsgDataDir, projName), sep='\t', header=True, index=True, float_format='%.2g')

# 10.5) Save the coordinates of the umap
rowIdx     = adataSubGroup.obs.index.values.tolist()
umapCordDF = pd.DataFrame(adataSubGroup.obsm['X_umap'],index=rowIdx, columns=['UMAP1','UMAP2', 'UMAP3'])
umapCordDF.to_csv("{0}/07_{1}_subCellType_assigned_{2}_UMAP_Coordinates.txt" .format(dataDir, projName, cluster_key), sep='\t', header=True, index=True, index_label="CellID")

# 10.6) Save the normalized, log transformed, batch and cell cycle corrected data
CellTypeDF                                = adata.to_df()
CellTypeDF['OriginalLouvainCluster']      = adata.obs['louvain_r1']
CellTypeDF['OriginalCellType']            = adata.obs['cellType']
CellTypeDF['SubCellLouvainCluster']       = adataSubGroup.obs['louvain_r1']
CellTypeDF['SubCellType']                 = adataSubGroup.obs['subCellType']
CellTypeDF.to_csv("{0}/07_{1}_scran_normalized_counts_annotation.txt" .format(dataDir, projName, cluster_key), sep='\t', header=True, index=True, index_label="CellID")

# 10.7) Marker genes & cluster annotation
# Annotation of cluster r_0.5 with known marker genes
markerDir = "{0}/SubCategories/markerDir".format(plotsDir); create_dir(markerDir)
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

  fig = plt.figure(figsize=(20,6*nrows))
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

  fig = plt.figure(figsize=(20,6*nrows))
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

# Get pairwise data dir
onesgDataDir  = "{0}/rankedGenes/SubCategories_oneVsRest".format(dataDir); create_dir(onesgDataDir)
namesDF    = pd.DataFrame(adataSubGroup.uns['rank_genes_{0}'.format(cluster_key)]['names'])
namesDF.to_csv("{0}/04_SubCategories_{1}_GeneNames.txt".format(onesgDataDir, projName), sep='\t', header=True, index=False, float_format='%.2g')
logfcDF    = pd.DataFrame(adataSubGroup.uns['rank_genes_{0}'.format(cluster_key)]['logfoldchanges'])
logfcDF.to_csv("{0}/04_SubCategories_{1}_logfoldchanges.txt".format(onesgDataDir, projName), sep='\t', header=True, index=False, float_format='%.2g')
pvalsadjDF = pd.DataFrame(adataSubGroup.uns['rank_genes_{0}'.format(cluster_key)]['pvals_adj'])
pvalsadjDF.to_csv("{0}/04_SubCategories_{1}_padj.txt".format(onesgDataDir, projName), sep='\t', header=True, index=False, float_format='%.2g')

# 10.7) Save the subCellType assigned adataSubGroup into a file
# Write the adataSubGroup and cadataSubGroup object to file
adataSubGroupfile  = "{0}/07_subCellType_assigned_{1}_adataSubGroup.h5ad" .format(dataDir, projName); adataSubGroup.write(adataSubGroupfile)
# # Read back the corrected adataSubGroup object
# adataSubGroupfile  = "{0}/07_subCellType_assigned_{1}_adataSubGroup.h5ad" .format(dataDir, projName); markerSubGroupadata = sc.read_h5ad(adataSubGroupfile)
# adataSubGroup = markerSubGroupadata.copy()

# Finished on 2020-05May-14 03:33
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
