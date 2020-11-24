# pwd
cd '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq'

# # Copy the H5 files
# for d in 3306_tumorA  3306_tumorB  3520_Antrum 3520_Corpus; 
# do 
#  cp -rv /media/rad/HDD2/temp_manec/02StomachBaseline/output/featureCounts/${d}/outs/filtered_feature_bc_matrix.h5 /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/manec/stomachBaseline/${d}_filtered_feature_bc_matrix.h5; 
# done;


##################################################################
ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)
##################################################################

# Loading the python libraries environment
%load scripts/load_python_modules.py

%load scripts/scrnaseq_module.py

# For filtering criteria
# https://krishnaswamylab.github.io/tutorial/load_filter_transform/
# For the below dataset, I would remove all cells with more than 25,000 UMI / cell 
# in fear they might represent doublets of cells. I will generally also remove all 
# cells with fewer than 500 reads per cell.

# System variables and directories
projName        = "stomachBaseline_mitGFP_3520_Antrum"
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/{0}".format(projName); create_dir("{0}".format(output_dir))
ccGenes_macosko = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/macosko_cell_cycle_genes_mmu.txt"
ccGenes_regev   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/regev_lab_cell_cycle_genes_mmu.txt"
minGenesPerCell = 50
minCountPerCell = 100
maxCountPerCell = 25000
minCellsPergene = 5
mtGenesFilter   = 0.35
rbGenesFilter   = 0.35
bname           = projName
plotsDir        = "{0}/plots".format(output_dir); create_dir(plotsDir)
dataDir         = "{0}/data".format(output_dir); create_dir(dataDir)
markerDir       = "{0}/markerDir".format(plotsDir); create_dir(markerDir)

# Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

# Get the logging parameters
sc.logging.print_memory_usage()
sc.logging.print_version_and_date()
sc.logging.print_header()
# Memory usage: current 1.01 GB, difference +1.01 GB
# Running Scanpy 1.6.0, on 2020-10-31 21:04.
# scanpy==1.6.0 anndata==0.7.4 umap==0.4.4 numpy==1.18.5 scipy==1.4.1 pandas==1.1.2 scikit-learn==0.23.1 statsmodels==0.11.1 python-igraph==0.8.2 louvain==0.7.0

# 1) Reading and performing QC on individual datasets
# 1.1) Reading the data in the anndata object individually
adata = sc.read_10x_h5('input/manec/stomachBaseline/3520_Antrum_filtered_feature_bc_matrix.h5')

# 1.2) Make the variable names unique and calculate some general qc-stats for genes and cells
adata.var_names_make_unique()

# Convert the sparse count matrices to dense represntation
adata.X = adata.X.toarray()

# 1.6) Calculate and plot QC covariates/metrices
rawadata = perform_qc(adata, plotsDir, bname, batch_key=None)
# - Shape (10000, 21828)
# - Plot unfiltered QC data
# Total number of cells: 10000
# Number of cells after min count filter ( 100  ): 10000
# Number of cells after max count filter (25000 ): 9956
# Number of cells after MT filter   (0.350000): 2809
# Number of cells after Ribo filter (0.350000): 2804
# Trying to set attribute `.obs` of view, copying.
# Number of cells after gene filter (  50  ): 2803
# Total number of genes: 21828
# Number of genes after minCellsPergene filter (  5   ): 13298
# - Filtered rawqcadata shape: (2803, 13298)

# Number of highly variable genes: 4924

# 1.7) Save the filtered raw adata into a file
# Write the adata object to file
adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata.write(adatafile)
# # Read back the filtered raw adata object
# adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata  = sc.read_h5ad(adatafile)

########################
# 2) Normalization
########################
scranadata = scran_norm(rawadata, plotsDir, bname)

# ########################
# # 3) Cell cycle correction
# ########################
cell_cycle_correction(scranadata, plotsDir, bname)

# 3.2) Save the normalized adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/02_norm_scran_ccc_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/02_norm_scran_ccc_adata.h5ad".format(dataDir, projName); scranadata = sc.read_h5ad(adatafile)
# adata = scranadata.copy()

########################
# 4) Clustering
########################
# 4.1) Perform clustering - using highly variable genes
sc.tl.leiden(adata                , key_added='leiden'     , random_state=2105)
sc.tl.leiden(adata, resolution=1  , key_added='leiden_r1'  , random_state=2105)
sc.tl.leiden(adata, resolution=1.5, key_added='leiden_r1.5', random_state=2105)
sc.tl.leiden(adata, resolution=2.0, key_added='leiden_r2'  , random_state=2105)

for i in np.linspace(0.1,0.9,9):
    try:
        sc.tl.leiden(adata, resolution=i, key_added='leiden_r{0}'.format(i), random_state=2105)
        print(adata.obs['leiden_r{0:0.1f}'.format(i)].value_counts())
    except:
        print("- Error in r: {0}".format(i))
sc.tl.leiden(adata, resolution=0.3, key_added='leiden_r0.3', random_state=2105)
sc.tl.leiden(adata, resolution=0.7, key_added='leiden_r0.7', random_state=2105)

# 4.2) Visualize the clustering and how this is reflected by different technical covariates
# UMAP
sc.pl.umap(adata, color=['leiden', 'leiden_r0.1', 'leiden_r0.2', 'leiden_r0.3', 'leiden_r0.4', 'leiden_r0.5', 'leiden_r0.6', 'leiden_r0.7', 'leiden_r0.8', 'leiden_r0.9', 'leiden_r1', 'leiden_r1.5', 'leiden_r2'], palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_{1}_clustering_all_leiden_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adata, color=['leiden', 'leiden_r0.1', 'leiden_r0.2', 'leiden_r0.3', 'leiden_r0.4', 'leiden_r0.5', 'leiden_r0.6', 'leiden_r0.7', 'leiden_r0.8', 'leiden_r0.9', 'leiden_r1', 'leiden_r1.5', 'leiden_r2'], palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/03_{1}_clustering_all_leiden_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
# TSNE
sc.pl.tsne(adata, color=['leiden', 'leiden_r0.1', 'leiden_r0.2', 'leiden_r0.3', 'leiden_r0.4', 'leiden_r0.5', 'leiden_r0.6', 'leiden_r0.7', 'leiden_r0.8', 'leiden_r0.9', 'leiden_r1', 'leiden_r1.5', 'leiden_r2'], palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_{1}_clustering_all_leiden_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 4.3) Select the resolution for clusters
leidenadata      = adata.copy()
cluster_key      = "leiden_r1"
cluster_bname    = "leiden_r1"
fig              = plt.figure(figsize=(26,10))
nElementsCluster = len(adata.obs[cluster_key].value_counts())
colPalette       = sc.pl.palettes.vega_20_scanpy
if 21 <= nElementsCluster <= 28:
  colPalette = sc.pl.palettes.vega_20
else:
  colPalette = sc.pl.palettes.godsnot_102

# 2D projection
ax = fig.add_subplot(2, 4, 1); sc.pl.tsne(adata, legend_loc='on data', ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(2, 4, 2);                  sc.pl.umap(adata   , legend_loc='on data', ax=ax, color=cluster_key   , palette=colPalette, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format(cluster_key))
ax = fig.add_subplot(2, 4, 3);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="log_counts", palette=colPalette, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="log_counts UMAP")
ax = fig.add_subplot(2, 4, 4);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="mt_frac"   , palette=colPalette, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="mt_frac UMAP")
# 3D projection
ax = fig.add_subplot(2, 4, 5); sc.pl.tsne(adata, ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(2, 4, 6, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color=cluster_key, palette=colPalette, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format(cluster_key))
ax = fig.add_subplot(2, 4, 7, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color="log_counts"   , palette=colPalette, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="log_counts UMAP")
ax = fig.add_subplot(2, 4, 8, projection='3d'); sc.pl.umap(adata  , legend_loc=None, ax=ax, color="mt_frac"   , palette=colPalette, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="mt_frac UMAP")
plt.tight_layout()
plt.savefig("{0}/03_{1}_{2}_sampleID_counts_mtfrac_UMAP.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 4.4) leiden UMAPs and TSNEs
fig = plt.figure(figsize=(24,14))
fig.suptitle("{0} UMAP".format(cluster_key))
ax = fig.add_subplot(2, 3, 1); sc.pl.tsne(adata, legend_loc='on data', ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.5, show=False)
ax = fig.add_subplot(2, 3, 2); sc.pl.umap(adata, legend_loc='on data', ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.5, show=False)
ax = fig.add_subplot(2, 3, 3, projection='3d'); sc.pl.umap(adata, ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
ax = fig.add_subplot(2, 3, 4); sc.pl.tsne(adata, legend_loc=None, ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.5, show=False)
ax = fig.add_subplot(2, 3, 5); sc.pl.umap(adata, ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.5, show=False)
plt.savefig("{0}/03_{1}_clustering_{2}_2D3D_UMAP_TSNE.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 4.5) Plot individual samples
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key=cluster_key, cluster_bname=cluster_bname, analysis_stage_num='03', analysis_stage='clustering', final_color_palette=colPalette)

# 4.6) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/03_LeidenClustering_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/03_LeidenClustering_{1}_adata.h5ad" .format(dataDir, projName); leidenadata  = sc.read_h5ad(adatafile)
# adata = leidenadata.copy()

########################
# 5) Marker genes identification & cluster annotation
########################
# 5.1.1) Calculate marker genes
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby=cluster_key, key_added='rank_genes_{0}'.format(cluster_key), n_genes=adata.shape[1])

# 5.1.2) Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_{0}'.format(cluster_key), fontsize=12, show=False)
plt.savefig("{0}/04_{1}_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# 5.2.1) Annotation of cluster r1 with known marker genes
subplot_title_fontsize = 12
subplot_title_width    = 50
cluster_key        = "leiden_r1"
cluster_bname      = "leiden_r1"

# 5.2.2) Get all the gene names in the adata object
genespresent = adata.var.index.values.tolist()

# 5.2.3) Read the marker genes into a pandas dataframe
marker_file            = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_V2.txt'
marker_genes_cellTypes = import_marker_genes_list(marker_file)
# Generate the UMAPs/TSNE for each marker categories
plot_manual_marker_list_genes(adata, markerDir, bname, cluster_key, genespresent, marker_genes_cellTypes, "stomach_marker_list_V2")
# Other marker gene visualization
plot_additional_marker_gene_visualization(adata, markerDir, bname, cluster_key, genespresent, marker_genes_cellTypes, "stomach_marker_list_V2", analysis_stage_num='04', analysis_stage='norm')

# 5.2.4) For mouse cell atlas marker genes
ma_marker_file         = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_mouse_cellatlas_V1.txt'
ma_marker_genes        = import_marker_genes_list(marker_file)
# Generate the UMAPs/TSNE for each marker categories
plot_manual_marker_list_genes(adata, markerDir, bname, cluster_key, genespresent, ma_marker_genes, "stomach_marker_list_mouse_cellatlas_V1")
# Other marker gene visualization
plot_additional_marker_gene_visualization(adata, markerDir, bname, cluster_key, genespresent, ma_marker_genes, "stomach_marker_list_mouse_cellatlas_V1", analysis_stage_num='04', analysis_stage='norm')

# 5.2.5) For custom sanger vector and GFP list
sangfp_marker_file     = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/sanger_gfp_marker_list_V1.txt'
sangfp_marker_genes        = import_marker_genes_list(sangfp_marker_file, 'GFP_hgMycIresCd2_V1')
# Generate the UMAPs/TSNE for each marker categories
plot_manual_marker_list_genes(adata, markerDir, bname, cluster_key, genespresent, sangfp_marker_genes, "GFP_hgMycIresCd2_V1")
# Other marker gene visualization
plot_additional_marker_gene_visualization(adata, markerDir, bname, cluster_key, genespresent, sangfp_marker_genes, "GFP_hgMycIresCd2_V1", analysis_stage_num='04', analysis_stage='norm')

# 5.3) Save the leiden information in external file
leidenDF = pd.DataFrame(adata.obs[cluster_key])
leidenDF.to_csv("{0}/03_{1}_leiden.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="cellId")

# 8.3) Dataframe of ranked genes
# Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
cluster_key_groups = adata.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = adata.obs[cluster_key].value_counts().to_dict()
rankGenesDir       = "{0}/rankedGenes/{1}".format(dataDir,cluster_bname); create_dir(rankGenesDir)
for g in cluster_key_groups:
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(adata.uns['rank_genes_{0}'.format(cluster_key)][n])[g]
  # Save dataframes
  ngDF.to_csv("{0}/03_{1}_rank_genes_{2}_{3}.txt".format(rankGenesDir, projName, cluster_bname, g), sep='\t', header=True, index=False, float_format='%.2g')

# 8.4) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); markeradata  = sc.read_h5ad(adatafile)
# adata = markeradata.copy()
# Finished on 2020-10-05 01:58:23
