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
input_h5_files  = 'input/manec/stomachBaseline/3520_Antrum_filtered_feature_bc_matrix.h5'
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/{0}".format(projName); create_dir("{0}".format(output_dir))
ccGenes_macosko = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/macosko_cell_cycle_genes_mmu.txt"
ccGenes_regev   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/regev_lab_cell_cycle_genes_mmu.txt"
marker_lists    = ["sanger_gfp_marker_list_V1.txt","01_manec_markers_various_celltypes_V1.txt", "03_manec_markers_endocrine_cells_V1.txt","04_manec_markers_progenitor_cells_V1.txt","02_manec_markers_nonpriority_cells_V1.txt","05_manec_markers_mousecellatlas_SL_JG_V1.txt"]
minGenesPerCell = 100
minCountPerCell = 200
maxCountPerCell = 25000
minCellsPergene = 10
mtGenesFilter   = 0.25
rbGenesFilter   = 0.25
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

# 1.1) Reading and performing QC on individual datasets
adata = import_datasets(input_h5_files, sampleIdDict=None, batch_key=None)

# 1.2) Calculate and plot QC covariates/metrices
rawadata = perform_qc(adata, plotsDir, bname, batch_key=None)

# 1.3) Save the filtered raw adata into a file
# Write the adata object to file
adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata.write(adatafile)
# # Read back the filtered raw adata object
# adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata  = sc.read_h5ad(adatafile)

########################
# 2) Normalization and Cell cycle correction
########################
# 2.1) Scran normalization
scranadata = scran_norm(rawadata, plotsDir, bname)

# 2.2) Cell cycle correction
cell_cycle_correction(scranadata, plotsDir, bname)
adata = scranadata.copy()

# 2.3) Save the normalized adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/02_norm_scran_ccc_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/02_norm_scran_ccc_adata.h5ad".format(dataDir, projName); scranadata = sc.read_h5ad(adatafile)
# adata = scranadata.copy()

########################
# 3) Clustering
########################
# 3.1) Calculate the umap and tsne with spreadout parameters
adata = calculate_umap_tsne(adata, num_neighbors=5, perplexity=10, early_exaggeration=5, random_state=2105)

# 3.2) Perform clustering - using highly variable genes
adata = calculate_plot_clustering(adata, plotsDir, bname, main_title = 'Clustering all resolution TSNE/UMAP', additional_features=None, analysis_stage_num='03', analysis_stage='clustering', color_palette=sc.pl.palettes.godsnot_102, clustering_algorithm='leiden', random_state = 2105)

# 3.3) Select the resolution for clusters
leidenadata      = adata.copy()
cluster_key      = "leiden_r1"
cluster_bname    = "leiden_r1"
nElementsCluster = len(adata.obs[cluster_key].value_counts())
colPalette       = get_color_palette(nElementsCluster)

# 3.4) Plot selected features
cluster_features = [cluster_key, 'log_counts', 'mt_frac', 'rb_frac']
plot_umap_tsne(adata, plotsDir, "{0}_{1}_counts_mtfrac_UMAP_TSNE".format(bname, cluster_bname), main_title = 'Leiden Counts mt/rb-frac TSNE/UMAP', features=cluster_features, analysis_stage_num='03', analysis_stage='clustering', color_palette=colPalette)

# 3.5) leiden UMAPs and TSNEs
plot_selected_cluster_umap_tsne(adata, plotsDir, bname, main_title = 'Leiden clustering', features=cluster_key, additional_features=None, analysis_stage_num='03', analysis_stage='clustering', color_palette=colPalette)

# 3.6) Plot individual samples
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key=cluster_key, cluster_bname=cluster_bname, analysis_stage_num='03', analysis_stage='clustering', final_color_palette=colPalette)

# 3.7) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/03_LeidenClustering_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/03_LeidenClustering_{1}_adata.h5ad" .format(dataDir, projName); leidenadata  = sc.read_h5ad(adatafile)
# adata = leidenadata.copy()

########################
# 4) Marker genes identification & cluster annotation
########################
# 4.1.1) Calculate marker genes
# Select the resolution for clusters
cluster_key        = "leiden_r1"
cluster_bname      = "leiden_r1"

# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby=cluster_key, key_added='rank_genes_{0}'.format(cluster_key), n_genes=adata.shape[1])

# 4.1.2) Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_{0}'.format(cluster_key), fontsize=12, show=False)
plt.savefig("{0}/04_{1}_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# 4.2.1) Annotate using SCSA a cell type annotation for single-cell RNA-seq data
marker_file        = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/04_manec_markers_progenitor_cells_V1_scsa.txt'
marker_list_name   = "manec_markers_progenitor_cells_V1_scsa"
cellTypeColumnName = "cellType"
adata = annotate_with_scsa(adata, dataDir, cluster_key, marker_file, projName, marker_list_name, cellTypeColumnName, analysis_stage_num='05', analysis_stage='scsa_clusters')

# 4.2.2) Plot selected features
nElementsCluster = len(adata.obs['CellType'].value_counts())
colPalette       = get_color_palette(nElementsCluster)
cluster_features = ['CellType']
plot_umap_tsne(adata, plotsDir, "{0}_{1}_CellType_UMAP_TSNE".format(bname, marker_list_name), main_title = "{0} Celltype TSNE/UMAP".format(marker_list_name), features=cluster_features, analysis_stage_num='05', analysis_stage='scsa_clusters', color_palette=colPalette)

# 4.3.1) Get all the gene names in the adata object
genespresent = adata.var.index.values.tolist()

# 4.3.2) Read the marker genes into a pandas dataframe
for marfile in marker_lists:
  marker_list_name       = get_file_info(marfile)[1]
  marker_file            = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/{0}".format(marfile)
  marker_genes_cellTypes = import_marker_genes_list(marker_file)

  # Generate the UMAPs/TSNE for each marker categories
  plot_manual_marker_list_genes(adata, markerDir, bname, cluster_key, genespresent, marker_genes_cellTypes, marker_list_name)
  # Other marker gene visualization
  plot_additional_marker_gene_visualization(adata, markerDir, bname, cluster_key, genespresent, marker_genes_cellTypes, marker_list_name, analysis_stage_num='04', analysis_stage='marker_genes_analysis')

# 4.3.3) Save the leiden information in external file
leidenDF = pd.DataFrame(adata.obs[cluster_key])
leidenDF.to_csv("{0}/03_{1}_leiden.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="cellId")

# 4.3.4) Dataframe of ranked genes
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

# 4.4) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); markeradata  = sc.read_h5ad(adatafile)
# adata = markeradata.copy()
# Finished on 2020-10-05 01:58:23

########################################################################################################################################################
# GFP positive cells only
genelist = ['GFP']
gfpadata = adata[adata[:,genelist].X.sum(1) > 0]
gfpadata[gfpadata.to_df()['GFP'] >=4]

# Recaluclate the TSNE and UMAP parameters
gfpadata = calculate_umap_tsne(gfpadata, num_neighbors=15, perplexity=30, random_state=2105)

# Plot the UMAPs and TSNEs
qcfeatures = ['log_counts', 'mt_frac', 'rb_frac']
plot_umap_tsne(gfpadata, plotsDir, bname, main_title = 'GFP positive cells only QC features', features=qcfeatures, analysis_stage_num='05', analysis_stage='GFP_pos_coutns_mt_rb_frac', color_palette=sc.pl.palettes.vega_20_scanpy)

# 4.1) Perform clustering - using highly variable genes
sc.tl.leiden(gfpadata                , key_added='leiden'     , random_state=2105)
sc.tl.leiden(gfpadata, resolution=1  , key_added='leiden_r1'  , random_state=2105)
sc.tl.leiden(gfpadata, resolution=1.5, key_added='leiden_r1.5', random_state=2105)
sc.tl.leiden(gfpadata, resolution=2.0, key_added='leiden_r2'  , random_state=2105)

for i in np.linspace(0.1,0.9,9):
    try:
        sc.tl.leiden(gfpadata, resolution=i, key_added='leiden_r{0}'.format(i), random_state=2105)
        print(gfpadata.obs['leiden_r{0:0.1f}'.format(i)].value_counts())
    except:
        print("- Error in r: {0}".format(i))
sc.tl.leiden(gfpadata, resolution=0.3, key_added='leiden_r0.3', random_state=2105)
sc.tl.leiden(gfpadata, resolution=0.7, key_added='leiden_r0.7', random_state=2105)


# 4.2) Visualize the clustering and how this is reflected by different technical covariates
cluster_features = ['leiden', 'leiden_r0.1', 'leiden_r0.2', 'leiden_r0.3', 'leiden_r0.4', 'leiden_r0.5', 'leiden_r0.6', 'leiden_r0.7', 'leiden_r0.8', 'leiden_r0.9', 'leiden_r1', 'leiden_r1.5', 'leiden_r2']
plot_umap_tsne(gfpadata, plotsDir, "{0}_clustering_all_leiden_UMAP_TSNE".format(bname), main_title = 'Clustering all leiden TSNE/UMAP', features=cluster_features, analysis_stage_num='05', analysis_stage='GFP_pos', color_palette=sc.pl.palettes.godsnot_102)

# 4.3) Select the resolution for clusters
leidengfpadata      = gfpadata.copy()
cluster_key      = "leiden_r2"
cluster_bname    = "leiden_r2"
fig              = plt.figure(figsize=(26,10))
nElementsCluster = len(gfpadata.obs[cluster_key].value_counts())
colPalette       = sc.pl.palettes.vega_20_scanpy
if 21 <= nElementsCluster <= 28:
  colPalette = sc.pl.palettes.vega_20
else:
  colPalette = sc.pl.palettes.godsnot_102

# Plot selected features
cluster_features = [cluster_key, 'log_counts', 'mt_frac', 'rb_frac']
plot_umap_tsne(gfpadata, plotsDir, "{0}_{1}_sampleID_counts_mtfrac_UMAP_TSNE".format(bname, cluster_bname), main_title = 'Leiden sampleID Counts mt/rb-frac TSNE/UMAP', features=cluster_features, analysis_stage_num='05', analysis_stage='GFP_pos_clustering', color_palette=colPalette)

# 4.4) leiden UMAPs and TSNEs
fig = plt.figure(figsize=(24,14))
fig.suptitle("{0} UMAP".format(cluster_key))
ax = fig.add_subplot(2, 3, 1); sc.pl.tsne(gfpadata, legend_loc='on data', ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.5, show=False)
ax = fig.add_subplot(2, 3, 2); sc.pl.umap(gfpadata, legend_loc='on data', ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.5, show=False)
ax = fig.add_subplot(2, 3, 3, projection='3d'); sc.pl.umap(gfpadata, ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
ax = fig.add_subplot(2, 3, 4); sc.pl.tsne(gfpadata, legend_loc=None, ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.5, show=False)
ax = fig.add_subplot(2, 3, 5); sc.pl.umap(gfpadata, ax=ax, color=cluster_key, palette=colPalette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.5, show=False)
plt.savefig("{0}/05_{1}_GFP_pos_clustering_{2}_2D3D_UMAP_TSNE.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 4.5) Plot individual samples
plot_individual_cluster_umap(gfpadata, plotsDir, bname, cluster_key=cluster_key, cluster_bname=cluster_bname, analysis_stage_num='03', analysis_stage='clustering', final_color_palette=colPalette)

# 4.6) Save the cellType assigned gfpadata into a file
# Write the gfpadata and cgfpadata object to file
gfpadatafile  = "{0}/05_GFP_pos_LeidenClustering_{1}_gfpadata.h5ad" .format(dataDir, projName); gfpadata.write(gfpadatafile)
# # Read back the corrected gfpadata object
# gfpadatafile  = "{0}/05_GFP_pos_LeidenClustering_{1}_gfpadata.h5ad" .format(dataDir, projName); leidengfpadata  = sc.read_h5ad(gfpadatafile)
# gfpadata = leidengfpadata.copy()

#############################################################
# 5.1.1) Calculate marker genes
gfpmarkerDir = "{0}/GFPpos".format(markerDir); create_dir(gfpmarkerDir)
# Calculate marker genes
sc.tl.rank_genes_groups(gfpadata, groupby=cluster_key, key_added='rank_genes_{0}'.format(cluster_key), n_genes=gfpadata.shape[1])

# 5.1.2) Plot marker genes
sc.pl.rank_genes_groups(gfpadata, key='rank_genes_{0}'.format(cluster_key), fontsize=12, show=False)
plt.savefig("{0}/05_{1}_GFP_pos_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# 5.2.1) Annotation of cluster r1 with known marker genes
subplot_title_fontsize = 12
subplot_title_width    = 50
cluster_key        = "leiden_r1"
cluster_bname      = "leiden_r1"

# 5.2.2) Get all the gene names in the gfpadata object
genespresent = gfpadata.var.index.values.tolist()

# 5.2.3) Read the marker genes into a pandas dataframe
marker_file            = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_V2.txt'
marker_genes_cellTypes = import_marker_genes_list(marker_file)
# Generate the UMAPs/TSNE for each marker categories
plot_manual_marker_list_genes(gfpadata, gfpmarkerDir, bname, cluster_key, genespresent, marker_genes_cellTypes, "stomach_marker_list_V2")
# Other marker gene visualization
plot_additional_marker_gene_visualization(gfpadata, gfpmarkerDir, bname, cluster_key, genespresent, marker_genes_cellTypes, "stomach_marker_list_V2", analysis_stage_num='05', analysis_stage='GFP_pos_markers')

# 5.2.4) For mouse cell atlas marker genes
ma_marker_file         = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_mouse_cellatlas_V1.txt'
ma_marker_genes        = import_marker_genes_list(marker_file)
# Generate the UMAPs/TSNE for each marker categories
plot_manual_marker_list_genes(gfpadata, gfpmarkerDir, bname, cluster_key, genespresent, ma_marker_genes, "stomach_marker_list_mouse_cellatlas_V1")
# Other marker gene visualization
plot_additional_marker_gene_visualization(gfpadata, gfpmarkerDir, bname, cluster_key, genespresent, ma_marker_genes, "stomach_marker_list_mouse_cellatlas_V1", analysis_stage_num='05', analysis_stage='GFP_pos_markers')

# 5.2.5) For mouse cell atlas marker genes
ma_marker_file         = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_comprehensive_V1.txt'
ma_marker_genes        = import_marker_genes_list(marker_file)
# Generate the UMAPs/TSNE for each marker categories
plot_manual_marker_list_genes(gfpadata, gfpmarkerDir, bname, cluster_key, genespresent, ma_marker_genes, "stomach_marker_list_comprehensive_V1")
# Other marker gene visualization
plot_additional_marker_gene_visualization(gfpadata, gfpmarkerDir, bname, cluster_key, genespresent, ma_marker_genes, "stomach_marker_list_comprehensive_V1", analysis_stage_num='05', analysis_stage='GFP_pos_markers')

# 5.2.6) For custom sanger vector and GFP list
sangfp_marker_file     = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/sanger_gfp_marker_list_V1.txt'
sangfp_marker_genes        = import_marker_genes_list(sangfp_marker_file, 'GFP_hgMycIresCd2_V1')
# Generate the UMAPs/TSNE for each marker categories
plot_manual_marker_list_genes(gfpadata, gfpmarkerDir, bname, cluster_key, genespresent, sangfp_marker_genes, "GFP_hgMycIresCd2_V1")
# Other marker gene visualization
plot_additional_marker_gene_visualization(gfpadata, gfpmarkerDir, bname, cluster_key, genespresent, sangfp_marker_genes, "GFP_hgMycIresCd2_V1", analysis_stage_num='05', analysis_stage='GFP_pos_markers')

# 5.3) Save the leiden information in external file
leidenDF = pd.DataFrame(gfpadata.obs[cluster_key])
leidenDF.to_csv("{0}/04_{1}_GFP_pos_leiden.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="cellId")

# 8.3) Dataframe of ranked genes
# Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
cluster_key_groups = gfpadata.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = gfpadata.obs[cluster_key].value_counts().to_dict()
rankGenesDir       = "{0}/rankedGenes/{1}".format(dataDir,cluster_bname); create_dir(rankGenesDir)
for g in cluster_key_groups:
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(gfpadata.uns['rank_genes_{0}'.format(cluster_key)][n])[g]
  # Save dataframes
  ngDF.to_csv("{0}/05_{1}_GFP_pos_rank_genes_{2}_{3}.txt".format(rankGenesDir, projName, cluster_bname, g), sep='\t', header=True, index=False, float_format='%.2g')

# 8.4) Save the cellType assigned gfpadata into a file
# Write the gfpadata and cgfpadata object to file
gfpadatafile  = "{0}/05_GFP_pos_markerGenes_{1}_gfpadata.h5ad" .format(dataDir, projName); gfpadata.write(gfpadatafile)
# # Read back the corrected gfpadata object
# gfpadatafile  = "{0}/05_GFP_pos_markerGenes_{1}_gfpadata.h5ad" .format(dataDir, projName); markergfpadata  = sc.read_h5ad(gfpadatafile)
# gfpadata = markergfpadata.copy()
# Finished on 2020-10-05 01:58:23

outputFileName =  "05_GFP_pos_markerGenes_{0}_GFP.txt" .format(projName); gfpadata.write(gfpadatafile)
save_adata_to_excel(gfpadata, dataDir, outputFileName, obs_additional_colnames=None, append_new_colnames=False, obs_colnames=None, subset_genes=['GFP'])