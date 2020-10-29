# pwd
cd '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq'

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
projName        = "allSamples_ohneB1079B1001_trvae"
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/{0}".format(projName); create_dir("{0}".format(output_dir))
ccGenes_macosko = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/macosko_cell_cycle_genes_mmu.txt"
ccGenes_regev   = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/annotations/regev_lab_cell_cycle_genes_mmu.txt"
minGenesPerCell = 50
minCountPerCell = 100
maxCountPerCell = 25000
minCellsPergene = 5
mtGenesFilter   = 0.25
rbGenesFilter   = 0.25
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
# 1.1.1) Reading the data in the anndata object individually
# Except: 'input/manec/tissue/bulk1079_mouse_filtered_feature_bc_matrix.h5',
#         'input/manec/tissue/bulk1001_mouse_filtered_feature_bc_matrix.h5'

# [AnnData object with n_obs × n_vars = 18273 × 55488
#  AnnData object with n_obs × n_vars =  1852 × 55488]
tissueFilenames = [
                    'input/manec/tissue/bulk997_mouse_filtered_feature_bc_matrix.h5',
                    'input/manec/tissue/stomach1001_mouse_filtered_feature_bc_matrix.h5',
                    'input/manec/tissue/bulk1018_mouse_filtered_feature_bc_matrix.h5'
                  ]
adatas          = [sc.read_10x_h5(f) for f in tissueFilenames]
adatas
#  AnnData object with n_obs × n_vars =  4334 × 55488
#  AnnData object with n_obs × n_vars =   904 × 55488
#  AnnData object with n_obs × n_vars =   897 × 55488]

# 1.1.2) Make variable names unique
for ad in adatas:
    print(ad.shape)
    ad.var_names_make_unique()

# 1.2) Get the dictionary of tissue ids
sampleIdDict =  {
                    '0':'bulk997', 
                    '1':'stomach1001',
                    '2':'bulk1018'
                }

# 1.3) Merge 10x datasets for different mices
adata = adatas[0].concatenate(adatas[1:])

# 1.4) Make variable names unique
adata.var_names_make_unique()

# Convert the sparse count matrices to dense represntation
adata.X = adata.X.toarray()

# 1.5) Add tissue id column for the batches
adata.obs['sampleID']  = adata.obs['batch'].map(sampleIdDict)
# In [16]: adata.obs
# Out[16]:
#                     batch sampleID
# AAACAAACAGCTATGA-0      0     S503
# AAACAAACCTACGAGC-0      0     S503

# 1.6) Calculate and plot QC covariates/metrices
rawadata = perform_qc(adata, plotsDir, bname)
# - Total number of cells                       : 6135
# - Number of cells after min count filter      : 6135
# - Number of cells after max count filter      : 6135
# - Number of cells after MT filter             : 5732
# - Number of cells after Ribo filter           : 4916
# - Number of cells after gene filter           : 4916

# - Total number of genes                       : 55488
# - Number of genes after minCellsPergene filter: 12923

# - Unfiltered rawqcadata shape                 : (6135, 55488)
# - Filtered rawqcadata shape                   : (4916, 12923)


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
# adata = scranadata.copy()

########################
# 4) Batch correction
########################
normadata = adata.copy()
rawadatas = adatas.copy()

# 5) Technical correction: Batch Correction using trVAE
# Source: https://github.com/theislab/trVAE/blob/master/examples/trVAE_Haber.ipynb

# 5.1) Normalizing & Extracting Top 2000 Highly Variable Genes
# One can use more genes but in order to train the network quickly, we will extract top 2000 genes. This can be done with normalize_hvg function in the tl module of trVAE package. The function accepts the following arguments:
# - adata: adata containing raw counts in its .X attribute.
# - target_sum: total counts per cell after normalization
# - size_factors: whether to normalize the adata and put total counts per cell in "size_factors" column of adata.obs (True is recommended).
# - scale_input: whether to scale the dataset after normalization (False is recommended).
# - logtrans_input: whether to log-transform the adata after normalization (True is recommended).
# - n_top_genes: number of highly variable genes to be selected after adata normalization.

condition_key = "sampleID"
adata2 = trvae.tl.normalize_hvg(adata, target_sum=1e4, size_factors=True, scale_input=False, logtrans_input=True, n_top_genes=adata.shape[1])
# adata2 = trvae.tl.normalize_hvg(adata, target_sum=1e4, size_factors=True, scale_input=False, logtrans_input=True, n_top_genes=4000)
adata2

# 5.2) Calculate number of batches
conditions = adata2.obs[condition_key].unique().tolist()

# 5.3) Create the network
# Some of network parameters:
# - x_dimension: size input features (necessary)
# - conditons: list of unique batches(studies) names
# - architecture: architecture of the network (optional)
# - output_activation: activation function of trVAE's last layer
# - alpha: coefficient of KL Divergence loss (optional)
# - beta: coefficient of MMD loss (optional)
# - eta: coefficient of reconstruction (MSE or SSE) loss (optional) can be one of the relu, leaky_relu, linear, ...
# - gene_names: list of gene names (adata.var_names.tolist())
# - loss_fn: trVAE's loss function (Has to be one of mse or sse)
network = trvae.models.trVAE(x_dimension=adata2.shape[1], architecture=[256,64], z_dimension=10, gene_names=adata2.var_names.tolist(), conditions=conditions, model_path=dataDir, alpha=0.0001, beta=50, eta=100, loss_fn='sse', output_activation='linear')

# 5.4) Training trVAE
# You can train scArches with train function with the following parameters:
# adata: Annotated dataset used for training and evaluating scArches.
# - condition_key: name of the column in obs matrix in adata which contains the batch_id for each sample.
# - n_epochs: number of epochs used to train scArches.
# - batch_size: number of sample used to sample as mini-batches in order to optimize scArches. Please NOTE that for MSE loss with MMD regularization batch sizes upper that 512 is highly recommended
# - save: whether to save scArches' model and configs after training phase or not.
# - retrain: if False and scArches' pretrained model exists in model_path, will restore scArches' weights. Otherwise will train and validate scArches on adata.
network.train(adata2, condition_key, train_size=0.75, n_epochs=500, batch_size=2048, early_stop_limit=300, lr_reducer=10, verbose=5, save=True)

# 5.5) Get corrected gene expression data
# we transfer all conditions to the batch labels with maximum number of samples. target_condition is the the condtion that you want your source adata be transformed to.
adata2.obs[condition_key].value_counts()
# bulk997        3462
# bulk1001        971
# stomach1001     781
# bulk1018        673
# Name: sampleID, dtype: int64
target_condition = adata2.obs[condition_key].value_counts().index[0] # 'bulk997'
corrected_adata = network.predict(adata2,condition_key,target_condition=target_condition)

# 5.6) Add the original var 
# corrected_adata.var = adata2.var.copy()
corrected_adata.var = adata2.var.copy()

# 5.7) Compute variable genes
# We first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.
sc.pp.highly_variable_genes(corrected_adata, flavor='cell_ranger')
print('\n','Number of highly variable genes: {:d}'.format(np.sum(corrected_adata.var['highly_variable'])))
# Number of highly variable genes: 3928

# 5.8) UMAP visualization of corrected gene expression
# adata = corrected_adata.copy()
sc.pp.pca(corrected_adata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(corrected_adata, random_state = 2105, n_neighbors=10)
sc.tl.umap(corrected_adata, random_state = 2105, n_components=3)
# sc.pp.neighbors(corrected_adata, n_neighbors=25)
# sc.tl.umap(corrected_adata, random_state = 2105, n_components=3)

fig = plt.figure(figsize=(16,13))
fig.suptitle('sampleID')
# 2D projection
ax = fig.add_subplot(2, 2, 1);                  sc.pl.umap(rawadata      , legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 2, 2);                  sc.pl.umap(corrected_adata, legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE UMAP")
# 3D projection 
ax = fig.add_subplot(2, 2, 3, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 2, 4, projection='3d'); sc.pl.umap(corrected_adata           , ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE UMAP")
plt.tight_layout()
plt.savefig("{0}/03_norm_TrVAE_batchCorrection_{1}_sampleID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Plot individual samples
plot_individual_cluster_umap(corrected_adata, plotsDir, bname, cluster_key='sampleID', cluster_bname='sampleID', analysis_stage_num='03', analysis_stage='scran_trvae')

# Plot tsne
sc.tl.tsne(corrected_adata   , random_state = 2105, n_pcs=50)
sc.tl.tsne(rawadata, random_state = 2105, n_pcs=50)

# SampleID 2D projection
fig = plt.figure(figsize=(16,8))
fig.suptitle('sampleID')
ax = fig.add_subplot(1, 2, 1); sc.pl.tsne(rawadata, legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw UMAP")
ax = fig.add_subplot(1, 2, 2); sc.pl.tsne(corrected_adata                    , ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE UMAP")
plt.tight_layout()
plt.savefig("{0}/03_norm_TrVAE_batchCorrection_{1}_sampleID_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 3.2) Save the normalized batch corrected adata into a file
adata = corrected_adata.copy()
# Write the adata and cadata object to file
adatafile  = "{0}/03_norm_TrVAE_batchCorrection_{1}_adata.h5ad" .format(dataDir, projName); corrected_adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/03_norm_TrVAE_batchCorrection_{1}_adata.h5ad" .format(dataDir, projName); trvaeBCadata  = sc.read_h5ad(adatafile)
# adata = trvaeBCadata.copy()

########################
# 5) Automated annotation using SingleR
########################
adataMatFile  = "{0}/04_{1}_scran_trvae_adata.matrix" .format(dataDir, projName);
adata.to_df().T.to_csv(adataMatFile, index=True, header=True, sep="\t")
os.system("Rscript scripts/R_annotate_cells_using_singleR.R -if={0} -of={0}".format(adataMatFile))
ImmGenAnnFile   = "{0}/04_{1}_scran_trvae_adata_SingleR_ImmGenRef.txt"     .format(dataDir, projName);
MouseRnaAnnFile = "{0}/04_{1}_scran_trvae_adata_SingleR_MouseRNAseqRef.txt".format(dataDir, projName);

# 5.1) Read in the SingleR annotated cell data
ImmGenAnnDF = pd.read_csv(ImmGenAnnFile, sep='\t', usecols=['cellIDs', 'labels'], index_col=['cellIDs'], quoting=3);
adata.obs['ImmGenLabels'] = ImmGenAnnDF['labels'].astype('category').values
MouseRnaseqAnnDF = pd.read_csv(MouseRnaAnnFile, sep='\t', usecols=['cellIDs', 'labels'], index_col=['cellIDs'], quoting=3);
adata.obs['MouseRnaseqLabels'] = MouseRnaseqAnnDF['labels'].astype('category').values

# 5.2) UMAPs
fig = plt.figure(figsize=(40,16)); c = 3
# 2D projection
ax = fig.add_subplot(2, c, 1);                  sc.pl.umap(adata, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="SampleID")
ax = fig.add_subplot(2, c, 2);                  sc.pl.umap(adata, ax=ax, color="ImmGenLabels"     , legend_loc='right margin', palette=sc.pl.palettes.godsnot_102, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="ImmGenLabels UMAP")
ax = fig.add_subplot(2, c, 3);                  sc.pl.umap(adata, ax=ax, color="MouseRnaseqLabels", legend_loc='right margin', palette=sc.pl.palettes.godsnot_102, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="MouseRnaseqLabels UMAP")
# 3D projection
ax = fig.add_subplot(2, c, 4, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="SampleID")
ax = fig.add_subplot(2, c, 5, projection='3d'); sc.pl.umap(adata                 , ax=ax, color="ImmGenLabels"     , palette=sc.pl.palettes.godsnot_102, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="ImmGenLabels UMAP")
ax = fig.add_subplot(2, c, 6, projection='3d'); sc.pl.umap(adata                 , ax=ax, color="MouseRnaseqLabels", palette=sc.pl.palettes.godsnot_102, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="MouseRnaseqLabels UMAP")
# Save plot
plt.tight_layout()
plt.savefig("{0}/04_{1}_sampleID_ImmGenLabels_MouseRnaseqLabels_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 5.3) Plot individual samples
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key='MouseRnaseqLabels', cluster_bname='MouseRnaseqLabels', analysis_stage_num='04', analysis_stage='singleR_UMAP', color_palette="godsnot_102")
# Save test marker gene plot
sc.pl.umap(adata, color=['MouseRnaseqLabels', 'Foxp1', 'Cd3d'], use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.tight_layout()
plt.savefig("{0}/04_{1}_MarkerGene_testList_MouseRnaseqLabels_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 5.4) TSNE SampleID 2D projection
fig = plt.figure(figsize=(16,8))
fig.suptitle('sampleID')
ax = fig.add_subplot(2, 2, 1); sc.pl.tsne(rawadata, legend_loc=None, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw TSNE")
ax = fig.add_subplot(2, 2, 2); sc.pl.tsne(adata                    , ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE TSNE")
ax = fig.add_subplot(2, 2, 3); sc.pl.tsne(adata                    , ax=ax, color="ImmGenLabels"     , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="ImmGenLabels TSNE")
ax = fig.add_subplot(2, 2, 4); sc.pl.tsne(adata                    , ax=ax, color="MouseRnaseqLabels", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="MouseRnaseqLabels TSNE")
plt.tight_layout()
plt.savefig("{0}/04_{1}_sampleID_ImmGenLabels_MouseRnaseqLabels_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 5.5) Save the singleR annotated adata into a file
# Write the adata object to file
adatafile  = "{0}/04_{1}_singleR_annotated_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/04_{1}_singleR_annotated_adata.h5ad" .format(dataDir, projName); singleRadata  = sc.read_h5ad(adatafile)
# adata = singleRadata.copy()

########################
# 6) TUMOR Analysis
########################
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

# 7.2) Plot UMAP and TSNE for all the tumor cells
adata.obs['AllTumor_hgMycIresCd2'] = humaniresmyc
fig = plt.figure(figsize=(24,6))
fig.suptitle('AllTumor_hgMycIresCd2')
# 2D projection
ax = fig.add_subplot(1, 3, 1);                  
sc.pl.umap(adata, legend_loc=None, ax=ax, color="AllTumor_hgMycIresCd2", color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(1, 3, 2, projection='3d'); 
sc.pl.umap(adata, ax=ax, color="AllTumor_hgMycIresCd2", color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, projection='3d', show=False)
# TSNE
ax = fig.add_subplot(1, 3, 3);
sc.pl.tsne(adata, legend_loc=None, ax=ax, color="AllTumor_hgMycIresCd2", color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, show=False)
plt.savefig("{0}/05_normTrVAE_{1}_AllTumor_hgMycIresCd2_CellIDs_UMAP_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

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
plt.savefig("{0}/05_normTrVAE_{1}_Tumor_hgMycIresCd2_individual_CellIDs_upsetPlot.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Plot all AllTumor_hgMycIresCd2 vector components in separate UMAP
fig = plt.figure(figsize=(16,6))
fig.suptitle('AllTumor_hgMycIresCd2 vector components')
# 2D projection
sc.pl.umap(adata, color=['AllTumor_hgMycIresCd2','humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2'], color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, show=False)
plt.savefig("{0}/05_normTrVAE_{1}_Tumor_hgMycIresCd2_Components_UMAP_2D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
# 3D projection
sc.pl.umap(adata, color=['AllTumor_hgMycIresCd2','humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2'], color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/05_normTrVAE_{1}_Tumor_hgMycIresCd2_Components_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# TSNE
fig = plt.figure(figsize=(16,6))
fig.suptitle('AllTumor_hgMycIresCd2 vector components')
# 2D projection
sc.pl.tsne(adata, color=['AllTumor_hgMycIresCd2','humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2'], color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.8, hspace=0.35, wspace=0.3, show=False)
plt.savefig("{0}/05_normTrVAE_{1}_Tumor_hgMycIresCd2_Components_TSNE_2D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# Get tumor cells per tissue (and in percentage)
tumorDF = adata.obs[['batch','humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2','AllTumor_hgMycIresCd2']].copy()
# Replace batch number with batch names
tumorDF.replace({'batch': sampleIdDict}, inplace=True)
# Remove index for groupby
tumorDF.reset_index(drop=True, inplace=True)

# Get the number of cells for each cluster in every tissue
ncellstumorDF = tumorDF.groupby(['batch']).sum()
ncellstumorDF.to_csv("{0}/03_{1}_tumorComponents_per_tissueid.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="tissueID")

########################
# 7) Clustering
########################
# 7.1) Perform clustering - using highly variable genes
sc.tl.louvain(adata                , key_added='louvain'     , random_state=2105)
sc.tl.louvain(adata, resolution=1  , key_added='louvain_r1'  , random_state=2105)
sc.tl.louvain(adata, resolution=1.5, key_added='louvain_r1.5', random_state=2105)
sc.tl.louvain(adata, resolution=2.0, key_added='louvain_r2'  , random_state=2105)

for i in np.linspace(0.1,0.9,9):
    try:
        sc.tl.louvain(adata, resolution=i, key_added='louvain_r{0}'.format(i), random_state=2105)
        print(adata.obs['louvain_r{0:0.1f}'.format(i)].value_counts())
    except:
        print("- Error in r: {0}".format(i))
sc.tl.louvain(adata, resolution=0.3, key_added='louvain_r0.3', random_state=2105)
sc.tl.louvain(adata, resolution=0.7, key_added='louvain_r0.7', random_state=2105)

# 7.2) Visualize the clustering and how this is reflected by different technical covariates
# UMAP
sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/05_{1}_clustering_all_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/05_{1}_clustering_all_louvain_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# TSNE
sc.pl.tsne(adata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/05_{1}_clustering_all_louvain_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

louvainadata = adata.copy()

cluster_key   = "louvain_r1"
cluster_bname = "louvain_r1"
fig = plt.figure(figsize=(32,8))
# 2D projection
ax = fig.add_subplot(2, 5, 1);                  sc.pl.umap(rawadata,                  ax=ax, color="sampleID"  , palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw tissueID UMAP")
ax = fig.add_subplot(2, 5, 2);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="sampleID"  , palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE tissueID UMAP")
ax = fig.add_subplot(2, 5, 3);                  sc.pl.umap(adata   , legend_loc='on data', ax=ax, color=cluster_key   , palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format(cluster_key))
ax = fig.add_subplot(2, 5, 4);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="log_counts", palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="log_counts UMAP")
ax = fig.add_subplot(2, 5, 5);                  sc.pl.umap(adata   , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="mt_frac UMAP")
# 3D projection
ax = fig.add_subplot(2, 5, 6, projection='3d'); sc.pl.umap(rawadata, legend_loc=None,  ax=ax, color="sampleID", palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw tissueID UMAP")
ax = fig.add_subplot(2, 5, 7, projection='3d'); sc.pl.umap(adata   , legend_loc=None,  ax=ax, color="sampleID", palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE tissueID UMAP")
ax = fig.add_subplot(2, 5, 8, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color=cluster_key, palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format(cluster_key))
ax = fig.add_subplot(2, 5, 9, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color="log_counts"   , palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="log_counts UMAP")
ax = fig.add_subplot(2, 5, 10, projection='3d'); sc.pl.umap(adata  , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="mt_frac UMAP")
plt.tight_layout()
plt.savefig("{0}/05_normTrVAE_{1}_{2}_sampleID_counts_mtfrac_UMAP.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Louvain UMAPs
fig = plt.figure(figsize=(16,6))
fig.suptitle("{0} UMAP".format(cluster_key))
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  
sc.pl.umap(adata, legend_loc='on data', ax=ax, color=cluster_key, palette=sc.pl.palettes.zeileis_28, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(1, 2, 2, projection='3d'); 
sc.pl.umap(adata, ax=ax, color=cluster_key, palette=sc.pl.palettes.zeileis_28, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/05_normTrVAE_{1}_clustering_{2}_UMAP_2D3D.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 5.3) Plot individual samples
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key=cluster_key, cluster_bname=cluster_bname, analysis_stage_num='05', analysis_stage='clustering')

# Louvain TSNEs
fig = plt.figure(figsize=(16,6))
fig.suptitle("{0} TSNE".format(cluster_key))
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  
sc.pl.tsne(adata, legend_loc='on data', ax=ax, color=cluster_key, palette=sc.pl.palettes.zeileis_28, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(1, 2, 2);                  
sc.pl.tsne(adata, ax=ax, color=cluster_key, palette=sc.pl.palettes.zeileis_28, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
plt.savefig("{0}/05_normTrVAE_{1}_clustering_{2}_TSNE_2D.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 7.2) Plot separate bar plots, coloured in by cluster annotation, for each tissue
# Convert palette into colormap
clcmap = ListedColormap(sc.pl.palettes.zeileis_28)
# Get the DF of tissue and clusters
clusterBatchDF = adata.obs[['batch','{0}'.format(cluster_key)]].copy()
# Replace batch number with batch names
clusterBatchDF.replace({'batch': sampleIdDict}, inplace=True)
# Remove index for groupby
clusterBatchDF.reset_index(drop=True, inplace=True)
# Get the number of cells for each cluster in every tissue
ncellsClusterBatchDF = clusterBatchDF.groupby(['batch','{0}'.format(cluster_key)]).size()
# Get the percent of cells for each cluster in every tissue 
pcellsClusterBatchDF = pd.crosstab(index=clusterBatchDF['batch'], columns=clusterBatchDF['{0}'.format(cluster_key)], values=clusterBatchDF['{0}'.format(cluster_key)], aggfunc='count', normalize='index')
# Plot the barplots
fig = plt.figure(figsize=(32,24)); fig.suptitle("Cells for each {0} in each tissue".format(cluster_key))
# plot numbers of cells
ax = fig.add_subplot(2, 2, 1); ncellsClusterBatchDF.unstack().plot(kind='barh', stacked=True, colormap=clcmap, ax=ax, legend=None, title="Number of cells")
# plot percent of cells
ax = fig.add_subplot(2, 2, 2); pcellsClusterBatchDF.plot(kind='barh',stacked=True, colormap=clcmap, ax=ax, title="% of cells")
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='{0}'.format(cluster_key), title_fontsize=12)

# Get the number of cells for each tissue in every cluster
nbatchPerClusterIdDF = clusterBatchDF.groupby(['{0}'.format(cluster_key),'batch']).size()
# Get the percent of cells for each tissue in every cluster 
pbatchPerClusterIdDF = pd.crosstab(index=clusterBatchDF['{0}'.format(cluster_key)], columns=clusterBatchDF['batch'], values=clusterBatchDF['batch'], aggfunc='count', normalize='index')
# Plot the barplots
ax = fig.add_subplot(2, 2, 3); nbatchPerClusterIdDF.unstack().plot(kind='barh', stacked=True, colormap=clcmap, ax=ax, legend=None, title="number of cells for each tissue in every cluster")
# plot percent of cells
ax = fig.add_subplot(2, 2, 4); pbatchPerClusterIdDF.plot(kind='barh',stacked=True, colormap=clcmap, ax=ax, title="% of cells for each tissue in every cluster")
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='{0}'.format(cluster_key), title_fontsize=12)

# Save plots in a 2x2 grid style
plt.tight_layout() # For non-overlaping subplots
plt.savefig("{0}/05_normTrVAE_{1}_clustering_{2}_tissueID_cluster_barplot.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')


########################
# 8) Marker genes identification & cluster annotation
########################
# 8.1.1) Calculate marker genes
# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby=cluster_key, key_added='rank_genes_{0}'.format(cluster_key), n_genes=adata.shape[1])

# 8.1.2) Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_{0}'.format(cluster_key), fontsize=12, show=False)
plt.savefig("{0}/06_normTrVAE_{1}_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# 8.2.1) Annotation of cluster r_0.5 with known marker genes
markerDir = "{0}/markerDir".format(plotsDir); create_dir(markerDir)
subplot_title_fontsize = 12
subplot_title_width    = 50

# # 8.2.2) Get all the gene names in the adata object
genespresent = adata.var.index.values.tolist()

# # 8.2.3) Read the marker genes into a pandas dataframe
marker_file  = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_V2.txt'
markersDF    = pd.read_csv(marker_file, sep="\t")
marker_genes = markersDF.groupby('CellLines')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()
marker_genes_cellTypes = markersDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()
# Generate the UMAPs/TSNE for each marker categories
plot_manual_marker_list_genes(adata, markerDir, bname, cluster_key, genespresent, marker_genes_cellTypes, "stomach_marker_list_V2")

# # 8.2.4) For mouse cell atlas marker genes
ma_marker_file       = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_mouse_cellatlas_V1.txt'
ma_markersDF         = pd.read_csv(ma_marker_file, sep="\t", header=0, index_col=None)
ma_marker_genes      = ma_markersDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()
# Generate the UMAPs/TSNE for each marker categories
plot_manual_marker_list_genes(adata, markerDir, bname, cluster_key, genespresent, ma_marker_genes, "stomach_marker_list_mouse_cellatlas_V1")

# # 8.3) Save the louvain information in external file
louvainsDF = pd.DataFrame(adata.obs[cluster_key])
louvainsDF.to_csv("{0}/03_{1}_louvains.txt".format(dataDir, projName), sep='\t', header=True, index=True, index_label="cellId")

# 8.2) Other marker gene visualization
# 8.2.1) Get clean marker dict
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

marker_list_name = "stomach_V2"
# 8.2.2) Dot plots: The dotplot visualization provides a compact way of showing per group, the fraction of cells expressing a gene (dot size) and the mean expression of the gene in those cell (color scale).
# The use of the dotplot is only meaningful when the counts matrix contains zeros representing no gene counts. dotplot visualization does not work for scaled or corrected matrices in which cero counts had been replaced by other values.
sc.pl.dotplot(adata, marker_genes_filtered_dict, groupby=cluster_key, log=True, figsize=(40,12), show=False, dendrogram=True)
plt.savefig("{0}/05_norm_{1}_{2}_31_marker_genes_{3}_dotplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.3) Matrix plots: The matrixplot shows the mean expression of a gene in a group by category as a heatmap. In contrast to dotplot, the matrix plot can be used with corrected and/or scaled counts. By default raw counts are used.
sc.pl.matrixplot(adata, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False,cmap='Reds',  figsize=(40,12), standard_scale='group', show=False)
plt.savefig("{0}/05_norm_{1}_{2}_31_marker_genes_{3}_scaled_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.matrixplot(adata, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False, cmap='Reds', figsize=(40,12), standard_scale='group', vmin=0.5, show=False)
plt.savefig("{0}/05_norm_{1}_{2}_31_marker_genes_{3}_scaled_vmin0_05_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.4) Tracksplots: The track plot shows the same information as the heatmap, but, instead of a color scale, the gene expression is represented by height.
ad = adata.copy()
# ad.raw.X.data = np.exp(ad.raw.X.data)
ax = sc.pl.tracksplot(ad, marker_genes_filtered_dict, groupby=cluster_key, log=True, dendrogram=True, show=False, figsize=(50,30))
plt.savefig("{0}/05_norm_{1}_{2}_31_marker_genes_{3}_tracksplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')

marker_list_name = "mouse_cellatlas"
# 8.2.5) Dot plots
sc.pl.dotplot(adata, ma_marker_genes_filtered_dict, groupby=cluster_key, log=True, figsize=(40,12), show=False, dendrogram=True)
plt.savefig("{0}/05_norm_{1}_{2}_32_marker_genes_{3}_dotplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.6) Matrix plots
sc.pl.matrixplot(adata, ma_marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False,cmap='Reds',  figsize=(40,12), standard_scale='group', show=False)
plt.savefig("{0}/05_norm_{1}_{2}_32_marker_genes_{3}_scaled_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.matrixplot(adata, ma_marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False, cmap='Reds', figsize=(40,12), standard_scale='group', vmin=0.5, show=False)
plt.savefig("{0}/05_norm_{1}_{2}_32_marker_genes_{3}_scaled_vmin0_05_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.7) Tracksplots
ax = sc.pl.tracksplot(ad, ma_marker_genes_filtered_dict, groupby=cluster_key, log=True, dendrogram=True, show=False, figsize=(50,30))
plt.savefig("{0}/05_norm_{1}_{2}_32_marker_genes_{3}_tracksplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')

# 8.3) Dataframe of ranked genes
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

# 8.4) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/04_markerGenes_{1}_adata.h5ad" .format(dataDir, projName); markeradata  = sc.read_h5ad(adatafile)
# adata = markeradata.copy()
# Finished on 2020-10-05 01:58:23

###############################################
# 7.2) Plot separate bar plots, coloured in by cluster annotation, for each tissue
# Get the DF of tissue and clusters
clusterTumorDF = adata.obs[['batch','{0}'.format(cluster_key),'AllTumor_hgMycIresCd2']].copy()
# Remove index for groupby
clusterTumorDF.reset_index(drop=True, inplace=True)
# Get the number of cells for each cluster in every tissue
# ncellsClusterTumorDF = clusterTumorDF.groupby(['{0}'.format(cluster_key),'AllTumor_hgMycIresCd2']).size()
# # Get the percent of cells for each cluster in every tissue 
# pcellsClusterTumorDF = ncellsClusterTumorDF.groupby(level=0).apply(lambda x: 100 * x / float(x.sum()))

# Get the percent of cells for each cluster in every tissue 
ncellsClusterTumorDF = clusterTumorDF.groupby(['{0}'.format(cluster_key),'AllTumor_hgMycIresCd2']).size().to_frame()
ncellsClusterTumorDF.rename(columns = {0:'cellCount'}, inplace=True)
ncellsClusterTumorDF['tumorPC'] = 100 * ncellsClusterTumorDF['cellCount'] / ncellsClusterTumorDF.groupby(cluster_key)['cellCount'].transform('sum')
# ## -- End pasted text --
# In [18]: ncellsClusterTumorDF
# Out[18]:
#                                   cellCount    tumorPC
# louvain_r1 AllTumor_hgMycIresCd2
# 0          0                            800  94.899170
#            1                             43   5.100830

# 1          0                            550  96.660808
#            1                             19   3.339192

# 2          0                            533  94.671403
#            1                             30   5.328597

# 3          0                            411  84.394251
#            1                             76  15.605749

# 4          0                            318  93.529412
#            1                             22   6.470588

# 5          0                            145  47.540984
#            1                            160  52.459016

# 6          0                            196  66.440678
#            1                             99  33.559322

# 7          0                            229  95.815900
#            1                             10   4.184100

# 8          0                            196  95.609756
#            1                              9   4.390244

# 9          0                             64  31.372549
#            1                            140  68.627451
           
# 10         0                            184  98.395722
#            1                              3   1.604278

# 11         0                            162  91.525424
#            1                             15   8.474576

# 12         0                            133  83.125000
#            1                             27  16.875000

# 13         0                            107  94.690265
#            1                              6   5.309735

# 14         0                              2   2.061856
#            1                             95  97.938144

# 15         0                              6   7.594937
#            1                             73  92.405063

# 16         0                             50  94.339623
#            1                              3   5.660377

# sc.pl.umap(adata, color=[cluster_key, 'Myc','Cd2','humanCd2'], use_raw=False, color_map=mymap, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, ncols=2, legend_loc='on data')

################################################################################################################
# 9) Only Tumor cells analysis
# Get the tumor adata
tumoradata = adata[adata.obs['AllTumor_hgMycIresCd2'] == 1].copy()

# Calculations for UMAP and TSNE
# UMAP
sc.pp.pca(tumoradata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(tumoradata, random_state = 2105)
# sc.pp.neighbors(tumoradata, random_state = 2105, n_neighbors=10)
sc.tl.umap(tumoradata, random_state = 2105, n_components=3)

# TSNE
sc.tl.tsne(tumoradata   , random_state = 2105, n_pcs=50)


# 7.1) Perform clustering - using highly variable genes
sc.tl.louvain(tumoradata                , key_added='louvain'     , random_state=2105)
sc.tl.louvain(tumoradata, resolution=1  , key_added='louvain_r1'  , random_state=2105)
sc.tl.louvain(tumoradata, resolution=1.5, key_added='louvain_r1.5', random_state=2105)
sc.tl.louvain(tumoradata, resolution=2.0, key_added='louvain_r2'  , random_state=2105)

for i in np.linspace(0.1,0.9,9):
    try:
        sc.tl.louvain(tumoradata, resolution=i, key_added='louvain_r{0}'.format(i), random_state=2105)
        print(tumoradata.obs['louvain_r{0:0.1f}'.format(i)].value_counts())
    except:
        print("- Error in r: {0}".format(i))
sc.tl.louvain(tumoradata, resolution=0.3, key_added='louvain_r0.3', random_state=2105)
sc.tl.louvain(tumoradata, resolution=0.7, key_added='louvain_r0.7', random_state=2105)

# 7.2) Visualize the clustering and how this is reflected by different technical covariates
# UMAP
sc.pl.umap(tumoradata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/07_{1}_tomor_clustering_all_louvain_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.umap(tumoradata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, projection='3d', show=False)
plt.savefig("{0}/07_{1}_tumor_clustering_all_louvain_UMAP_3D.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

# TSNE
sc.pl.tsne(tumoradata, color=['louvain', 'louvain_r0.1', 'louvain_r0.2', 'louvain_r0.3', 'louvain_r0.4', 'louvain_r0.5', 'louvain_r0.6', 'louvain_r0.7', 'louvain_r0.8', 'louvain_r0.9', 'louvain_r1', 'louvain_r1.5', 'louvain_r2'], palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.savefig("{0}/07_{1}_tumor_clustering_all_louvain_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')


cluster_key   = "louvain_r1"
cluster_bname = "louvain_r1"
fig = plt.figure(figsize=(38,12))
# 2D projection
ax = fig.add_subplot(2, 5, 1);                  sc.pl.umap(rawadata,                  ax=ax, color="sampleID"  , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw tissueID UMAP")
ax = fig.add_subplot(2, 5, 2);                  sc.pl.umap(tumoradata   , legend_loc=None, ax=ax, color="sampleID"  , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE tissueID UMAP")
ax = fig.add_subplot(2, 5, 3);                  sc.pl.umap(tumoradata   , legend_loc='on data', ax=ax, color=cluster_key   , palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format(cluster_key))
ax = fig.add_subplot(2, 5, 4);                  sc.pl.umap(tumoradata   , legend_loc=None, ax=ax, color="log_counts", palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="log_counts UMAP")
ax = fig.add_subplot(2, 5, 5);                  sc.pl.umap(tumoradata   , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="mt_frac UMAP")
# 3D projection
ax = fig.add_subplot(2, 5, 6, projection='3d'); sc.pl.umap(tumoradata,    legend_loc=None,  ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw tissueID UMAP")
ax = fig.add_subplot(2, 5, 7, projection='3d'); sc.pl.umap(tumoradata   , legend_loc=None,  ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE tissueID UMAP")
ax = fig.add_subplot(2, 5, 8, projection='3d'); sc.pl.umap(tumoradata   ,                   ax=ax, color=cluster_key, palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format(cluster_key))
ax = fig.add_subplot(2, 5, 9, projection='3d'); sc.pl.umap(tumoradata   , legend_loc=None, ax=ax, color="log_counts"   , palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="log_counts UMAP")
ax = fig.add_subplot(2, 5, 10, projection='3d'); sc.pl.umap(tumoradata  , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.zeileis_28, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="mt_frac UMAP")
plt.tight_layout()
plt.savefig("{0}/07_{1}_tumor_{2}_sampleID_counts_mtfrac_UMAP.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Louvain UMAPs
fig = plt.figure(figsize=(16,6))
fig.suptitle("{0} UMAP".format(cluster_key))
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  
sc.pl.umap(tumoradata, legend_loc='on data', ax=ax, color=cluster_key, palette=sc.pl.palettes.zeileis_28, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
# 3D projection
ax = fig.add_subplot(1, 2, 2, projection='3d'); 
sc.pl.umap(tumoradata, ax=ax, color=cluster_key, palette=sc.pl.palettes.zeileis_28, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
plt.savefig("{0}/07_{1}_tumor_clustering_{2}_UMAP_2D3D.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 5.3) Plot individual samples
plot_individual_cluster_umap(tumoradata, plotsDir, bname, cluster_key=cluster_key, cluster_bname=cluster_bname, analysis_stage_num='07', analysis_stage='tumor_clustering')

# Louvain TSNEs
fig = plt.figure(figsize=(16,6))
fig.suptitle("{0} TSNE".format(cluster_key))
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  
sc.pl.tsne(tumoradata, legend_loc='on data', ax=ax, color=cluster_key, palette=sc.pl.palettes.zeileis_28, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
ax = fig.add_subplot(1, 2, 2);                  
sc.pl.tsne(tumoradata, ax=ax, color=cluster_key, palette=sc.pl.palettes.zeileis_28, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
plt.savefig("{0}/07_{1}_tumor_clustering_{2}_TSNE_2D.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 7.2) Plot separate bar plots, coloured in by cluster annotation, for each tissue
# Convert palette into colormap
clcmap = ListedColormap(sc.pl.palettes.zeileis_28)
# Get the DF of tissue and clusters
clusterBatchDF = tumoradata.obs[['batch','{0}'.format(cluster_key)]].copy()
# Replace batch number with batch names
clusterBatchDF.replace({'batch': sampleIdDict}, inplace=True)
# Remove index for groupby
clusterBatchDF.reset_index(drop=True, inplace=True)
# Get the number of cells for each cluster in every tissue
ncellsClusterBatchDF = clusterBatchDF.groupby(['batch','{0}'.format(cluster_key)]).size()
# Get the percent of cells for each cluster in every tissue 
pcellsClusterBatchDF = pd.crosstab(index=clusterBatchDF['batch'], columns=clusterBatchDF['{0}'.format(cluster_key)], values=clusterBatchDF['{0}'.format(cluster_key)], aggfunc='count', normalize='index')
# Plot the barplots
fig = plt.figure(figsize=(32,24)); fig.suptitle("Cells for each {0} in each tissue".format(cluster_key))
# plot numbers of cells
ax = fig.add_subplot(2, 2, 1); ncellsClusterBatchDF.unstack().plot(kind='barh', stacked=True, colormap=clcmap, ax=ax, legend=None, title="Number of cells")
# plot percent of cells
ax = fig.add_subplot(2, 2, 2); pcellsClusterBatchDF.plot(kind='barh',stacked=True, colormap=clcmap, ax=ax, title="% of cells")
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='{0}'.format(cluster_key), title_fontsize=12)

# Get the number of cells for each tissue in every cluster
nbatchPerClusterIdDF = clusterBatchDF.groupby(['{0}'.format(cluster_key),'batch']).size()
# Get the percent of cells for each tissue in every cluster 
pbatchPerClusterIdDF = pd.crosstab(index=clusterBatchDF['{0}'.format(cluster_key)], columns=clusterBatchDF['batch'], values=clusterBatchDF['batch'], aggfunc='count', normalize='index')
# Plot the barplots
ax = fig.add_subplot(2, 2, 3); nbatchPerClusterIdDF.unstack().plot(kind='barh', stacked=True, colormap=clcmap, ax=ax, legend=None, title="number of cells for each tissue in every cluster")
# plot percent of cells
ax = fig.add_subplot(2, 2, 4); pbatchPerClusterIdDF.plot(kind='barh',stacked=True, colormap=clcmap, ax=ax, title="% of cells for each tissue in every cluster")
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='{0}'.format(cluster_key), title_fontsize=12)

# Save plots in a 2x2 grid style
plt.tight_layout() # For non-overlaping subplots
plt.savefig("{0}/07_{1}_tumor_clustering_{2}_tissueID_cluster_barplot.png".format(plotsDir, bname, cluster_bname) , bbox_inches='tight', dpi=175); plt.close('all')

# 8.1.1) Calculate marker genes
# Calculate marker genes
sc.tl.rank_genes_groups(tumoradata, groupby=cluster_key, key_added='rank_genes_{0}'.format(cluster_key), n_genes=adata.shape[1])

# 8.1.2) Plot marker genes
sc.pl.rank_genes_groups(tumoradata, key='rank_genes_{0}'.format(cluster_key), fontsize=12, show=False)
plt.savefig("{0}/07_{1}_{2}_tumor_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# # 8.2.2) Get all the gene names in the adata object
tumorgenespresent = tumoradata.var.index.values.tolist()

# Generate the UMAPs/TSNE for each marker categories
tumormarkerDir = "{0}/markerDir/tumor".format(plotsDir); create_dir(tumormarkerDir)
subplot_title_fontsize = 12
subplot_title_width    = 50
plot_manual_marker_list_genes(tumoradata, tumormarkerDir, bname, cluster_key, tumorgenespresent, marker_genes_cellTypes, "stomach_marker_list_V2")
plot_manual_marker_list_genes(tumoradata, tumormarkerDir, bname, cluster_key, tumorgenespresent, ma_marker_genes, "stomach_marker_list_mouse_cellatlas_V1")

# 8.2.1) Get clean marker dict
adata_expressed_genes = tumoradata.var.index.tolist()
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

# Get the dendrogram reflect the same cluster_key
sc.tl.dendrogram(tumoradata, groupby=cluster_key)

marker_list_name = "stomach_V2"
# 8.2.2) Dot plots: The dotplot visualization provides a compact way of showing per group, the fraction of cells expressing a gene (dot size) and the mean expression of the gene in those cell (color scale).
# The use of the dotplot is only meaningful when the counts matrix contains zeros representing no gene counts. dotplot visualization does not work for scaled or corrected matrices in which cero counts had been replaced by other values.
sc.pl.dotplot(tumoradata, marker_genes_filtered_dict, groupby=cluster_key, log=True, figsize=(40,12), show=False, dendrogram=True)
plt.savefig("{0}/007_{1}_{2}_tumor_marker_genes_{3}_dotplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.3) Matrix plots: The matrixplot shows the mean expression of a gene in a group by category as a heatmap. In contrast to dotplot, the matrix plot can be used with corrected and/or scaled counts. By default raw counts are used.
sc.pl.matrixplot(tumoradata, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False,cmap='Reds',  figsize=(40,12), standard_scale='group', show=False)
plt.savefig("{0}/07_{1}_{2}_tumor_marker_genes_{3}_scaled_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.matrixplot(tumoradata, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False, cmap='Reds', figsize=(40,12), standard_scale='group', vmin=0.5, show=False)
plt.savefig("{0}/07_{1}_{2}_tumor_marker_genes_{3}_scaled_vmin0_05_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.4) Tracksplots: The track plot shows the same information as the heatmap, but, instead of a color scale, the gene expression is represented by height.
ad = tumoradata.copy()
# ad.raw.X.data = np.exp(ad.raw.X.data)
ax = sc.pl.tracksplot(ad, marker_genes_filtered_dict, groupby=cluster_key, log=True, dendrogram=True, show=False, figsize=(50,30))
plt.savefig("{0}/07_{1}_{2}_tumor_marker_genes_{3}_tracksplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')

marker_list_name = "mouse_cellatlas"
# 8.2.5) Dot plots
sc.pl.dotplot(tumoradata, ma_marker_genes_filtered_dict, groupby=cluster_key, log=True, figsize=(40,12), show=False, dendrogram=True)
plt.savefig("{0}/07_{1}_{2}_tumor_marker_genes_{3}_dotplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.6) Matrix plots
sc.pl.matrixplot(tumoradata, ma_marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False,cmap='Reds',  figsize=(40,12), standard_scale='group', show=False)
plt.savefig("{0}/07_{1}_{2}_tumor_marker_genes_{3}_scaled_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
sc.pl.matrixplot(tumoradata, ma_marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False, cmap='Reds', figsize=(40,12), standard_scale='group', vmin=0.5, show=False)
plt.savefig("{0}/07_{1}_{2}_tumor_marker_genes_{3}_scaled_vmin0_05_matrixplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
# 8.2.7) Tracksplots
ax = sc.pl.tracksplot(ad, ma_marker_genes_filtered_dict, groupby=cluster_key, log=True, dendrogram=True, show=False, figsize=(50,30))
plt.savefig("{0}/07_{1}_{2}_tumor_marker_genes_{3}_tracksplot.png".format(plotsDir, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')

# 8.3) Dataframe of ranked genes
# Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
cluster_key        = "louvain_r1.5"
cluster_bname      = "louvain_r1_5"
cluster_key_groups = tumoradata.obs[cluster_key].cat.categories.tolist()
cluster_cell_count = tumoradata.obs[cluster_key].value_counts().to_dict()
rankGenesDir       = "{0}/rankedGenes/{1}/tumor".format(dataDir,cluster_bname); create_dir(rankGenesDir)
for g in cluster_key_groups:
  ngDF = pd.DataFrame()
  for n in ['names', 'scores', 'logfoldchanges',  'pvals', 'pvals_adj']:
    ngDF[n] = pd.DataFrame(tumoradata.uns['rank_genes_{0}'.format(cluster_key)][n])[g]
  # Save dataframes
  ngDF.to_csv("{0}/04_{1}_tumor_rank_genes_{2}_{3}.txt".format(rankGenesDir, projName, cluster_bname, g), sep='\t', header=True, index=False, float_format='%.2g')

# 8.4) Save the cellType assigned tumoradata into a file
# Write the tumoradata and ctumoradata object to file
tumoradatafile  = "{0}/05_tumor_markerGenes_{1}_tumoradata.h5ad" .format(dataDir, projName); tumoradata.write(tumoradatafile)
# # Read back the corrected tumoradata object
# tumoradatafile  = "{0}/05_tumor_markerGenes_{1}_tumoradata.h5ad".format(dataDir, projName); markertumoradata  = sc.read_h5ad(tumoradatafile)
# tumoradata = markertumoradata.copy()
# Finished on 2020-10-05 01:58:23