# pwd
cd '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq'

##################################################################
ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)
##################################################################
# https://scanpy-tutorials.readthedocs.io/en/latest/spatial/integration-scanorama.html

# Loading the python libraries environment
%load scripts/load_python_modules.py


%load scripts/scrnaseq_module.py

# For filtering criteria
# https://krishnaswamylab.github.io/tutorial/load_filter_transform/
# For the below dataset, I would remove all cells with more than 25,000 UMI / cell 
# in fear they might represent doublets of cells. I will generally also remove all 
# cells with fewer than 500 reads per cell.

# System variables and directories
projName        = "tregCNS"
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/{0}/allSamples_ohneS515_scanorama".format(projName); create_dir("{0}".format(output_dir))
ccGenes_macosko = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/macosko_cell_cycle_genes_mmu.txt"
ccGenes_regev   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/regev_lab_cell_cycle_genes_mmu.txt"
minGenesPerCell = 200
minCountPerCell = 200
maxCountPerCell = 20000 
minCellsPergene = 25
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
# Memory usage: current 1.04 GB, difference +1.04 GB
# Running Scanpy 1.4.6, on 2020-07-19 23:32.
# scanpy==1.4.6 anndata==0.7.3 umap==0.4.4 numpy==1.19.0rc2 scipy==1.4.1 pandas==0.25.3 scikit-learn==0.23.1 statsmodels==0.11.1 python-igraph==0.8.2 louvain==0.7.0

# 1) Reading and performing QC on individual datasets
# 1.1) Reading the data in the anndata object individually
# adata = sc.read_10x_mtx('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S509')
# S503/  S504/  S505/  S508/  S509/  S511/  S512/  S514/  S515/  S516/  S517/  S518/  S519/
SampleFileNames = [
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S503',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S504',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S505',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S508',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S509',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S511',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S512',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S514',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S516',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S517',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S518',
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S519'
                  ]
adatas          = [sc.read_10x_mtx(f) for f in SampleFileNames]
adatas

# 1.2) Get the dictionary of tissue ids
sampleIDDict  =  {'0' :'S503', '1' :'S504', '2' :'S505', '3' :'S508', '4' :'S509', '5' :'S511', '6' :'S512', '7' :'S514', '8' :'S516', '9' :'S517', '10':'S518', '11':'S519'}
conditionDict =  {'0' :'Naive', '1' :'Naive', '2' :'Naive', '3' :'Peak', '4' :'Peak', '5' :'Peak', '6' :'Peak', '7' :'Remission', '8' :'Remission', '9' :'Remission', '10':'Remission', '11':'Remission'}

# 1.3) Merge 10x datasets for different mices
adata = adatas[0].concatenate(adatas[1:])

# 1.4) Make variable names unique
adata.var_names_make_unique()

# Convert the sparse count matrices to dense represntation
adata.X = adata.X.toarray()

# 1.5) Add tissue id column for the batches
adata.obs['sampleID']  = adata.obs['batch'].map(sampleIDDict)
adata.obs['condition'] = adata.obs['batch'].map(conditionDict)
# In [16]: adata.obs
# Out[16]:
#                     batch sampleID
# AAACAAACAGCTATGA-0      0     S503
# AAACAAACCTACGAGC-0      0     S503

# 1.6) Calculate and plot QC covariates/metrices
rawadata = perform_qc(adata, plotsDir, bname)
# - Unfiltered rawqcadata shape: (28632, 55471)
#     - Total number of cells: 28632
#     - Number of cells after min count filter: 26139
#     - Number of cells after max count filter: 26116
#     - Number of cells after MT filter  : 26098
#     - Number of cells after Ribo filter: 25329
#     - Number of cells after gene filter: 24891

#     - Total number of genes: 55471
#     - Number of genes after minCellsPergene filter: 16617
#     - Number of highly variable genes: 2969
# - Filtered rawqcadata shape: (24891, 16617)

# Plot individual clusters
plot_individual_cluster_umap(rawadata, plotsDir, "{0}_filtered".format(bname), cluster_key='sampleID', cluster_bname='sampleID')

# 1.7) Save the filtered raw adata into a file
# Write the adata object to file
adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata.write(adatafile)
# # Read back the filtered raw adata object
# adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata  = sc.read_h5ad(adatafile)

########################
# 2) Normalization
########################
# cpmadata   = cpm_norm(rawadata, plotsDir, bname)
# scranadata = scran_norm(rawadata, plotsDir, bname)
# Perform a clustering for scran normalization in clusters
adata    = rawadata.copy() 
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
size_factors = computeSumFactors(data_mat, clusters=input_groups)
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
adata.raw  = adata
scranadata = adata.copy()
# ########################
# # 3) Cell cycle correction
# ########################
cell_cycle_correction(adata, plotsDir, bname)

# 3.2) Save the normalized adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/02_norm_scran_ccc_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/02_norm_scran_ccc_adata.h5ad".format(dataDir, projName); scranadata = sc.read_h5ad(adatafile)
# normadata = adata.copy()

########################
# 4) Batch correction
########################
normadata = adata.copy()
rawadatas = adatas.copy()

# 4.1) Technical correction: Batch Correction using scanorama
qcall_s = scanorama_bc(adata, plotsDir, bname, batchkey='batch')
# Add to the AnnData object
adata.obsm["SC"] = qcall_s

# corrected_data = adata.copy()
# UMAP visualization of corrected gene expression
# adata = corrected_data.copy()
# sc.pp.neighbors(corrected_data, n_neighbors=25)
# sc.tl.umap(corrected_data, random_state = 2105, n_components=3)

fig = plt.figure(figsize=(36,16)); c = 4
# 2D projection
ax = fig.add_subplot(2, c, 1);                  sc.pl.umap(rawadata, ax=ax, color="sampleID" , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw SampleID")
ax = fig.add_subplot(2, c, 2);                  sc.pl.umap(rawadata, ax=ax, color="condition", legend_loc='right margin', palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw Condition")
ax = fig.add_subplot(2, c, 3);                  sc.pl.umap(adata   , ax=ax, color="sampleID" , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Scanorama SampleID")
ax = fig.add_subplot(2, c, 4);                  sc.pl.umap(adata   , ax=ax, color="condition", legend_loc='right margin', palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Scanorama Condition")
# 3D projection
ax = fig.add_subplot(2, c, 5, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="sampleID"      , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw SampleID")
ax = fig.add_subplot(2, c, 6, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="condition"     , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw Condition")
ax = fig.add_subplot(2, c, 7, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Scanorama SampleID")
ax = fig.add_subplot(2, c, 8, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color="condition"     , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Scanorama Condition")
# Save plot
plt.tight_layout()
plt.savefig("{0}/03_norm_scanorama_batchCorrection_{1}_sampleID_condition_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Plot individual samples
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key='sampleID' , cluster_bname='sampleID' , analysis_stage_num='03', analysis_stage='scrannorm_scanorama_batchCorrection')
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key='condition', cluster_bname='condition', analysis_stage_num='03', analysis_stage='scrannorm_scanorama_batchCorrection')

# Plot tsne
sc.tl.tsne(adata   , random_state = 2105, perplexity=35, n_pcs=50)
sc.tl.tsne(rawadata, random_state = 2105, perplexity=35, n_pcs=50)

# sampleID
fig = plt.figure(figsize=(16,8))
fig.suptitle('sampleID')
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  sc.pl.tsne(rawadata, legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw TSNE")
ax = fig.add_subplot(1, 2, 2);                  sc.pl.tsne(adata                    , ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Scanorama TSNE")
plt.tight_layout()
plt.savefig("{0}/03_norm_scanorama_batchCorrection_{1}_sampleID_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Condition
fig = plt.figure(figsize=(16,8))
fig.suptitle('condition')
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  sc.pl.tsne(rawadata, legend_loc=None, ax=ax, color="condition", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw TSNE")
ax = fig.add_subplot(1, 2, 2);                  sc.pl.tsne(adata                    , ax=ax, color="condition", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Scanorama TSNE")
plt.tight_layout()
plt.savefig("{0}/03_norm_scanorama_batchCorrection_{1}_condition_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 3.2) Save the normalized batch corrected adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/03_norm_scanorama_batchCorrection_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/03_norm_scanorama_batchCorrection_{1}_adata.h5ad" .format(dataDir, projName); scanoramaBCadata  = sc.read_h5ad(adatafile)

########################
# 5) Automated annotation using SingleR
########################

adataMatFile  = "{0}/04_{1}_scran_scanorama_adata.matrix" .format(dataDir, projName);
adata.to_df().T.to_csv(adataMatFile, index=True, header=True, sep="\t")
os.system("Rscript scripts/R_annotate_cells_using_singleR.R -if={0} -of={0}".format(adataMatFile))
ImmGenAnnFile   = "{0}/04_{1}_scran_scanorama_adata_SingleR_ImmGenRef.txt"     .format(dataDir, projName);
MouseRnaAnnFile = "{0}/04_{1}_scran_scanorama_adata_SingleR_MouseRNAseqRef.txt".format(dataDir, projName);

########################
# 6) Clustering
########################
# 6.1) Read in the SingleR annotated cell data
ImmGenAnnDF = pd.read_csv(ImmGenAnnFile, sep='\t', usecols=['cellIDs', 'labels'], index_col=['cellIDs'], quoting=3);
adata.obs['ImmGenLabels'] = ImmGenAnnDF['labels'].astype('category').values
MouseRnaseqAnnDF = pd.read_csv(MouseRnaAnnFile, sep='\t', usecols=['cellIDs', 'labels'], index_col=['cellIDs'], quoting=3);
adata.obs['MouseRnaseqLabels'] = MouseRnaseqAnnDF['labels'].astype('category').values

# 6.2) UMAPs
fig = plt.figure(figsize=(64,16)); c = 4
# 2D projection
ax = fig.add_subplot(2, c, 1);                  sc.pl.umap(adata, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="SampleID")
ax = fig.add_subplot(2, c, 2);                  sc.pl.umap(adata, ax=ax, color="condition"        , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Condition")
ax = fig.add_subplot(2, c, 3);                  sc.pl.umap(adata, ax=ax, color="ImmGenLabels"     , palette=sc.pl.palettes.godsnot_102, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="ImmGenLabels UMAP")
ax = fig.add_subplot(2, c, 4);                  sc.pl.umap(adata, ax=ax, color="MouseRnaseqLabels", palette=sc.pl.palettes.godsnot_102, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="MouseRnaseqLabels UMAP")
# 3D projection
ax = fig.add_subplot(2, c, 5, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="SampleID")
ax = fig.add_subplot(2, c, 6, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="condition"         , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Condition")
ax = fig.add_subplot(2, c, 7, projection='3d'); sc.pl.umap(adata                 , ax=ax, color="ImmGenLabels"     , palette=sc.pl.palettes.godsnot_102, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="ImmGenLabels UMAP")
ax = fig.add_subplot(2, c, 8, projection='3d'); sc.pl.umap(adata                 , ax=ax, color="MouseRnaseqLabels", palette=sc.pl.palettes.godsnot_102, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="MouseRnaseqLabels UMAP")
# Save plot
plt.tight_layout()
plt.savefig("{0}/04_{1}_sampleID_condition__ImmGenLabels_MouseRnaseqLabels_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 6.3) Plot individual samples
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key='MouseRnaseqLabels', cluster_bname='MouseRnaseqLabels', analysis_stage_num='04', analysis_stage='singleR_UMAP', color_palette="godsnot_102")
# Save test marker gene plot
sc.pl.umap(adata, color=['MouseRnaseqLabels', 'Foxp1', 'Cd3d'], use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.tight_layout()
plt.savefig("{0}/04_{1}_MarkerGene_testList_MouseRnaseqLabels_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 6.4) Save the singleR annotated adata into a file
# Write the adata object to file
adatafile  = "{0}/04_{1}_singleR_annotated_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/03_norm_all_batchCorrection_adata.h5ad" .format(dataDir, projName); normadata  = sc.read_h5ad(adatafile)
# normadata = adata.copy()










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

# 4.3) Visualizations

# Calculations for the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata, random_state = 2105, use_rep = "SC")
sc.tl.umap(adata, random_state = 2105, n_components=3)

# Plot visualizations
# Visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/03_{1}_clustering_all_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/03_{1}_clustering_all_louvain_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
