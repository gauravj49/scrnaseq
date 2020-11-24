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
conditionDict =  {'0' :'Naive', '1' :'Naive', '2' :'Naive', '3' :'Peak', '4' :'Peak', '5' :'Peak', '6' :'Peak', '7' :'Remission', '8' :'Remission', '9' :'Remission', '10':'Remission', 'Remission':'Remission'}

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
# normadata = adata.copy()

########################
# 4) Batch correction
########################
normadata = adata.copy()
rawadatas = adatas.copy()

# 4.1) Technical correction: Batch Correction using scanorama
qcall_s = scanorama_bc(qcadata, plotsDir, bname, batchkey='sampleID')
# Add to the AnnData object
scanoramaadata = adata.copy()
scanoramaadata.obsm["SC"] = qcall_s

# 5.1) Normalizing & Extracting Top 2000 Highly Variable Genes
# One can use more genes but in order to train the network quickly, we will extract top 2000 genes. This can be done with normalize_hvg function in the tl module of trVAE package. The function accepts the following arguments:
# - adata: adata containing raw counts in its .X attribute.
# - target_sum: total counts per cell after normalization
# - size_factors: whether to normalize the adata and put total counts per cell in "size_factors" column of adata.obs (True is recommended).
# - scale_input: whether to scale the dataset after normalization (False is recommended).
# - logtrans_input: whether to log-transform the adata after normalization (True is recommended).
# - n_top_genes: number of highly variable genes to be selected after adata normalization.

condition_key = "sampleID"
adata2 = trvae.tl.normalize_hvg(adata, target_sum=1e4, size_factors=True, scale_input=False, logtrans_input=True, n_top_genes=2000)
adata2

# Calculate number of batches
conditions = adata2.obs[condition_key].unique().tolist()

# Create the network
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

# Training trVAE
# You can train scArches with train function with the following parameters:

# adata: Annotated dataset used for training and evaluating scArches.
# - condition_key: name of the column in obs matrix in adata which contains the batch_id for each sample.
# - n_epochs: number of epochs used to train scArches.
# - batch_size: number of sample used to sample as mini-batches in order to optimize scArches. Please NOTE that for MSE loss with MMD regularization batch sizes upper that 512 is highly recommended
# - save: whether to save scArches' model and configs after training phase or not.
# - retrain: if False and scArches' pretrained model exists in model_path, will restore scArches' weights. Otherwise will train and validate scArches on adata.
network.train(adata2, condition_key, train_size=0.8, n_epochs=500, batch_size=2048, early_stop_limit=300, lr_reducer=20, verbose=5, save=True, )

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
sc.pp.neighbors(adata, random_state = 2105, use_rep = "mmd_latent")
sc.tl.umap(adata, random_state = 2105, n_components=3)

# # Getting corrected latent adata
# # if you use trVAE for batch-removal we recommend to use z Latent space computed using get_latent function This function has the following parameters:
# # - adata: Annotated dataset to be transformed to latent space
# # - batch_key: Name of the column in obs matrix in adata which contains the study for each sample.
# latent_adata = network.get_latent(adata2, condition_key, return_z=True)
# latent_adata

# # Get corrected gene expression data¶
# # we transfer all conditions to the batch labels with maximum number of samples. target_condition is the the condtion that you want your source adata be transformed to
adata.obs[condition_key].value_counts()
target_condition = adata.obs[condition_key].value_counts().index[0]
corrected_data   = network.predict(adata,condition_key,target_condition=target_condition)

# # Compute variable genes
# # We first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.
# sc.pp.highly_variable_genes(corrected_data, flavor='cell_ranger')
# print('\n','Number of highly variable genes: {:d}'.format(np.sum(corrected_data.var['highly_variable'])))

# # UMAP visualization of corrected gene expression
# # adata = corrected_data.copy()
# sc.pp.neighbors(corrected_data, n_neighbors=25)
# sc.tl.umap(corrected_data, random_state = 2105, n_components=3)

fig = plt.figure(figsize=(16,13))
fig.suptitle('sampleID')
# 2D projection
ax = fig.add_subplot(2, 2, 1);                  sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 2, 2);                  sc.pl.umap(corrected_data   , legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE UMAP")
# 3D projection 
ax = fig.add_subplot(2, 2, 3, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 2, 4, projection='3d'); sc.pl.umap(corrected_data                    , ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE UMAP")
plt.tight_layout()
plt.savefig("{0}/03_norm_TrVAE_batchCorrection_{1}_sampleID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Plot individual samples
plot_individual_cluster_umap(corrected_data, plotsDir, bname, cluster_key='sampleID', cluster_bname='sampleID', analysis_stage_num='03', analysis_stage='scrannorm_scanorama_batchCorrection')


# Plot tsne
sc.tl.tsne(corrected_data, random_state = 2105, perplexity=35, n_pcs=50)
sc.tl.tsne(rawadata, random_state = 2105, perplexity=35, n_pcs=50)

fig = plt.figure(figsize=(16,8))
fig.suptitle('sampleID')
# 2D projection
ax = fig.add_subplot(1, 2, 1);                  sc.pl.tsne(rawadata, legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw UMAP")
ax = fig.add_subplot(1, 2, 2);                  sc.pl.tsne(corrected_data           , ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE UMAP")
plt.tight_layout()
plt.savefig("{0}/03_norm_TrVAE_batchCorrection_{1}_sampleID_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 3.2) Save the normalized batch corrected adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/03_norm_TrVAE_batchCorrection_{1}_adata.h5ad" .format(dataDir, projName); corrected_data.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/03_norm_TrVAE_batchCorrection_{1}_adata.h5ad" .format(dataDir, projName); trvaeBCadata  = sc.read_h5ad(adatafile)

########################
# 5) Automated annotation using SingleR
########################

adataMatFile  = "{0}/04_{1}_scran_trvae_adata.matrix" .format(dataDir, projName);
adata.to_df().T.to_csv(adataMatFile, index=True, header=True, sep="\t")
os.system("Rscript scripts/R_annotate_cells_using_singleR.R -if={0} -of={0}".format(adataMatFile))
ImmGenAnnFile   = "{0}/04_{1}_scran_trvae_adata_SingleR_ImmGenRef.txt"     .format(dataDir, projName);
MouseRnaAnnFile = "{0}/04_{1}_scran_trvae_adata_SingleR_MouseRNAseqRef.txt".format(dataDir, projName);

########################
# 6) Clustering
########################
# 6.1) Read in the SingleR annotated cell data
ImmGenAnnDF = pd.read_csv(ImmGenAnnFile, sep='\t', usecols=['cellIDs', 'labels'], index_col=['cellIDs'], quoting=3);
adata.obs['ImmGenLabels'] = ImmGenAnnDF['labels'].astype('category').values
MouseRnaseqAnnDF = pd.read_csv(MouseRnaAnnFile, sep='\t', usecols=['cellIDs', 'labels'], index_col=['cellIDs'], quoting=3);
adata.obs['MouseRnaseqLabels'] = MouseRnaseqAnnDF['labels'].astype('category').values

# 6.2) UMAPs
fig = plt.figure(figsize=(40,16)); c = 3
# 2D projection
ax = fig.add_subplot(2, c, 1);                  sc.pl.umap(adata, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="SampleID")
ax = fig.add_subplot(2, c, 2);                  sc.pl.umap(adata, ax=ax, color="ImmGenLabels"     , legend_loc='right margin', palette=p, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="ImmGenLabels UMAP")
ax = fig.add_subplot(2, c, 3);                  sc.pl.umap(adata, ax=ax, color="MouseRnaseqLabels", legend_loc='right margin', palette=sc.pl.palettes.godsnot_102, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="MouseRnaseqLabels UMAP")
# 3D projection
ax = fig.add_subplot(2, c, 4, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="SampleID")
ax = fig.add_subplot(2, c, 5, projection='3d'); sc.pl.umap(adata                 , ax=ax, color="ImmGenLabels"     , palette=sc.pl.palettes.godsnot_102, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="ImmGenLabels UMAP")
ax = fig.add_subplot(2, c, 6, projection='3d'); sc.pl.umap(adata                 , ax=ax, color="MouseRnaseqLabels", palette=sc.pl.palettes.godsnot_102, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="MouseRnaseqLabels UMAP")
# Save plot
plt.tight_layout()
plt.savefig("{0}/04_{1}_sampleID_ImmGenLabels_MouseRnaseqLabels_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

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
