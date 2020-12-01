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
projName        = "tregCNS"
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/{0}/allSamples_ohneS515_trVAE".format(projName); create_dir("{0}".format(output_dir))
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
tissueFilenames = [
                    '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S503',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S504',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S505',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S508',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S509',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S511',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S512',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S514',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S516',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S517',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S518',       '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/tregCNS/10Xcompatible/S519'
                  ]
adatas          = [sc.read_10x_mtx(f) for f in tissueFilenames]
adatas

# 1.2) Get the dictionary of tissue ids
sampleIdDict  =  {'0' :'S503', '1' :'S504', '2' :'S505', '3' :'S508', '4' :'S509', '5' :'S511', '6' :'S512', '7' :'S514', '8' :'S516', '9' :'S517', '10':'S518', '11':'S519'}
conditionDict =  {'0' :'Naive', '1' :'Naive', '2' :'Naive', '3' :'Peak', '4' :'Peak', '5' :'Peak', '6' :'Peak', '7' :'Remission', '8' :'Remission', '9' :'Remission', '10':'Remission', '11':'Remission'}

# 1.3) Merge 10x datasets for different mices
adata = adatas[0].concatenate(adatas[1:])

# 1.4) Make variable names unique
adata.var_names_make_unique()

# Convert the sparse count matrices to dense represntation
adata.X = adata.X.toarray()

# 1.5) Add tissue id column for the batches
adata.obs['sampleID']  = adata.obs['batch'].map(sampleIdDict)
adata.obs['condition'] = adata.obs['batch'].map(conditionDict)
# In [16]: adata.obs
# Out[16]:
#                     batch sampleID
# AAACAAACAGCTATGA-0      0     S503
# AAACAAACCTACGAGC-0      0     S503

# 1.6) Calculate and plot QC covariates/metrices
rawadata = perform_qc(adata, plotsDir, bname, 15, 25)
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
sc.pp.pca(corrected_adata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(corrected_adata, random_state = 2105, n_neighbors=10)
sc.tl.umap(corrected_adata, random_state = 2105, n_components=3)

fig = plt.figure(figsize=(16,13))
fig.suptitle('sampleID')
# 2D projection
ax = fig.add_subplot(2, 2, 1);                  sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 2, 2);                  sc.pl.umap(corrected_adata   , legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE UMAP")
# 3D projection 
ax = fig.add_subplot(2, 2, 3, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw UMAP")
ax = fig.add_subplot(2, 2, 4, projection='3d'); sc.pl.umap(corrected_adata                    , ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE UMAP")
plt.tight_layout()
plt.savefig("{0}/03_norm_TrVAE_batchCorrection_{1}_sampleID_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Plot individual samples
plot_individual_cluster_umap(corrected_adata, plotsDir, bname, cluster_key='sampleID', cluster_bname='sampleID', analysis_stage_num='03', analysis_stage='scrannorm_TrVAE_batchCorrection')
plot_individual_cluster_umap(corrected_adata, plotsDir, bname, cluster_key='condition', cluster_bname='condition', analysis_stage_num='03', analysis_stage='scrannorm_TrVAE_batchCorrection')

# Plot tsne
sc.tl.tsne(adata   , random_state = 2105, perplexity=25, n_pcs=50)
sc.tl.tsne(rawadata, random_state = 2105, perplexity=25, n_pcs=50)

# SampleID 2D projection
fig = plt.figure(figsize=(16,8))
fig.suptitle('sampleID')
ax = fig.add_subplot(1, 2, 1); sc.pl.tsne(rawadata, legend_loc=None, ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw UMAP")
ax = fig.add_subplot(1, 2, 2); sc.pl.tsne(adata                    , ax=ax, color="sampleID", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE UMAP")
plt.tight_layout()
plt.savefig("{0}/03_norm_TrVAE_batchCorrection_{1}_sampleID_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Condition 2D projection
fig = plt.figure(figsize=(16,8))
fig.suptitle('condition')
ax = fig.add_subplot(1, 2, 1); sc.pl.tsne(rawadata, legend_loc=None, ax=ax, color="condition", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw TSNE")
ax = fig.add_subplot(1, 2, 2); sc.pl.tsne(adata                    , ax=ax, color="condition", palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE TSNE")
plt.tight_layout()
plt.savefig("{0}/03_norm_TrVAE_batchCorrection_{1}_condition_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Plot sampleID and Conditions for raw and bcdata
fig = plt.figure(figsize=(36,16)); c = 4
# 2D projection
ax = fig.add_subplot(2, c, 1);                  sc.pl.umap(rawadata, ax=ax, color="sampleID" , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw SampleID")
ax = fig.add_subplot(2, c, 2);                  sc.pl.umap(rawadata, ax=ax, color="condition", legend_loc='right margin', palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Raw Condition")
ax = fig.add_subplot(2, c, 3);                  sc.pl.umap(adata   , ax=ax, color="sampleID" , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE SampleID")
ax = fig.add_subplot(2, c, 4);                  sc.pl.umap(adata   , ax=ax, color="condition", legend_loc='right margin', palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="TrVAE Condition")
# 3D projection
ax = fig.add_subplot(2, c, 5, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="sampleID"      , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw SampleID")
ax = fig.add_subplot(2, c, 6, projection='3d'); sc.pl.umap(rawadata, legend_loc=None, ax=ax, color="condition"     , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Raw Condition")
ax = fig.add_subplot(2, c, 7, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE SampleID")
ax = fig.add_subplot(2, c, 8, projection='3d'); sc.pl.umap(adata   , legend_loc=None, ax=ax, color="condition"     , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="TrVAE Condition")
# Save plot
plt.tight_layout()
plt.savefig("{0}/03_norm_trVAE_batchCorrection_{1}_sampleID_condition_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')


# 3.2) Save the normalized batch corrected adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/03_norm_TrVAE_batchCorrection_{1}_adata.h5ad" .format(dataDir, projName); corrected_adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/03_norm_TrVAE_batchCorrection_{1}_adata.h5ad" .format(dataDir, projName); trvaeBCadata  = sc.read_h5ad(adatafile)

########################
# 5) Automated annotation using SingleR
########################

# 5.1) Perform automated annotation using SingleR
adataMatFile  = "{0}/04_{1}_scran_trvae_adata.matrix" .format(dataDir, projName);
adata.to_df().T.to_csv(adataMatFile, index=True, header=True, sep="\t")
os.system("Rscript scripts/R_annotate_cells_using_singleR.R -if={0} -of={0}".format(adataMatFile))
ImmGenAnnFile   = "{0}/04_{1}_scran_trvae_adata_SingleR_ImmGenRef.txt"     .format(dataDir, projName);
MouseRnaAnnFile = "{0}/04_{1}_scran_trvae_adata_SingleR_MouseRNAseqRef.txt".format(dataDir, projName);

# 5.2) Read in the SingleR annotated cell data
ImmGenAnnDF = pd.read_csv(ImmGenAnnFile, sep='\t', usecols=['cellIDs', 'labels'], index_col=['cellIDs'], quoting=3);
adata.obs['ImmGenLabels'] = ImmGenAnnDF['labels'].astype('category').values
MouseRnaseqAnnDF = pd.read_csv(MouseRnaAnnFile, sep='\t', usecols=['cellIDs', 'labels'], index_col=['cellIDs'], quoting=3);
adata.obs['MouseRnaseqLabels'] = MouseRnaseqAnnDF['labels'].astype('category').values

# 5.3) Plot UMAPs
fig = plt.figure(figsize=(52,16)); c = 4
# 2D projection
ax = fig.add_subplot(2, c, 1);                  sc.pl.umap(adata, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="SampleID")
ax = fig.add_subplot(2, c, 2);                  sc.pl.umap(adata, ax=ax, color="condition"         , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Conditions")
ax = fig.add_subplot(2, c, 3);                  sc.pl.umap(adata, ax=ax, color="ImmGenLabels"     , legend_loc='right margin', palette=sc.pl.palettes.zeileis_28, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="ImmGenLabels UMAP")
ax = fig.add_subplot(2, c, 4);                  sc.pl.umap(adata, ax=ax, color="MouseRnaseqLabels", legend_loc='right margin', palette=sc.pl.palettes.godsnot_102, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="MouseRnaseqLabels UMAP")
# 3D projection
ax = fig.add_subplot(2, c, 5, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="SampleID")
ax = fig.add_subplot(2, c, 6, projection='3d'); sc.pl.umap(adata, legend_loc=None, ax=ax, color="condition"         , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="Conditions")
ax = fig.add_subplot(2, c, 7, projection='3d'); sc.pl.umap(adata                 , ax=ax, color="ImmGenLabels"     , palette=sc.pl.palettes.godsnot_102, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="ImmGenLabels UMAP")
ax = fig.add_subplot(2, c, 8, projection='3d'); sc.pl.umap(adata                 , ax=ax, color="MouseRnaseqLabels", palette=sc.pl.palettes.godsnot_102, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="MouseRnaseqLabels UMAP")
# Save plot
plt.tight_layout()
plt.savefig("{0}/04_{1}_sampleID_condition_ImmGenLabels_MouseRnaseqLabels_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# 5.4) Plot individual samples
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key='ImmGenLabels', cluster_bname='ImmGenLabels', analysis_stage_num='04', analysis_stage='singleR_UMAP', color_palette="zeileis_28")
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key='MouseRnaseqLabels', cluster_bname='MouseRnaseqLabels', analysis_stage_num='04', analysis_stage='singleR_UMAP', color_palette="godsnot_102")

# 5.5) Plot TSNEs
fig = plt.figure(figsize=(52,16)); c = 4
# 2D projection
ax = fig.add_subplot(2, c, 1);                  sc.pl.tsne(adata, ax=ax, color="sampleID"         , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="SampleID")
ax = fig.add_subplot(2, c, 2);                  sc.pl.tsne(adata, ax=ax, color="condition"         , palette=sc.pl.palettes.vega_20, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="Conditions")
ax = fig.add_subplot(2, c, 3);                  sc.pl.tsne(adata, ax=ax, color="ImmGenLabels"     , legend_loc='right margin', palette=sc.pl.palettes.zeileis_28, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="ImmGenLabels UMAP")
ax = fig.add_subplot(2, c, 4);                  sc.pl.tsne(adata, ax=ax, color="MouseRnaseqLabels", legend_loc='right margin', palette=sc.pl.palettes.godsnot_102, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="MouseRnaseqLabels UMAP")
# Save plot
plt.tight_layout()
plt.savefig("{0}/04_{1}_sampleID_condition_ImmGenLabels_MouseRnaseqLabels_TSNE.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# # Save test marker gene plot
# sc.pl.umap(adata, color=['MouseRnaseqLabels', 'Foxp1', 'Cd3d'], use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
# plt.tight_layout()
# plt.savefig("{0}/04_{1}_MarkerGene_testList_MouseRnaseqLabels_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')

# Save test marker gene plot
sc.pl.umap(adata, color=['louvain_r1', 'Foxp3', 'Cd3d', 'Cd4', 'Cd8a', 'Ptprc'], legend_loc='on data', use_raw=False, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, show=False)
plt.tight_layout()
plt.savefig("{0}/04_{1}_MarkerGene_testList2_louvain_r1_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=100); plt.close('all')


# 6.4) Save the singleR annotated adata into a file
# Write the adata object to file
adatafile  = "{0}/04_{1}_singleR_annotated_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # # Read back the corrected adata object
# adatafile  = "{0}/04_{1}_singleR_annotated_adata.h5ad" .format(dataDir, projName); singleradata  = sc.read_h5ad(adatafile)
# adata = singleradata.copy()


#######################
# 6) Clustering
########################
# 3.1) Calculate the umap and tsne with spreadout parameters
adata = calculate_umap_tsne(adata, num_neighbors=5, perplexity=10, early_exaggeration=5, random_state=2105)

# 3.2) Perform clustering - using highly variable genes
adata = calculate_plot_clustering(adata, plotsDir, bname, main_title = 'Clustering all resolution TSNE/UMAP', additional_features=None, analysis_stage_num='05', analysis_stage='clustering', color_palette=sc.pl.palettes.godsnot_102, clustering_algorithm='leiden', random_state = 2105)

# 3.3) Select the resolution for clusters
leidenadata      = adata.copy()
cluster_key      = "leiden_r1.5"
cluster_bname    = "leiden_r1.5"
nElementsCluster = len(adata.obs[cluster_key].value_counts())
colPalette       = get_color_palette(nElementsCluster)

# 3.4) Plot selected features
cluster_features = [cluster_key, 'sampleID', 'condition','log_counts', 'mt_frac', 'rb_frac']
plot_umap_tsne(adata, plotsDir, "{0}_{1}_counts_mtfrac_UMAP_TSNE".format(bname, cluster_bname), main_title = 'Leiden Counts mt/rb-frac TSNE/UMAP', features=cluster_features, analysis_stage_num='05', analysis_stage='clustering', color_palette=colPalette)

# 3.5) leiden UMAPs and TSNEs
plot_selected_cluster_umap_tsne(adata, plotsDir, bname, main_title = 'Leiden clustering', features=cluster_key, additional_features=['sampleID', 'condition'], analysis_stage_num='05', analysis_stage='clustering', color_palette=colPalette)

# 3.6) Plot individual samples
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key=cluster_key, cluster_bname=cluster_bname, analysis_stage_num='05', analysis_stage='clustering', final_color_palette=colPalette)

# 3.7) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/05_LeidenClustering_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/05_LeidenClustering_{1}_adata.h5ad" .format(dataDir, projName); leidenadata  = sc.read_h5ad(adatafile)
# adata = leidenadata.copy()

########################
# 4) Marker genes identification & cluster annotation
########################
# 4.1.1) Calculate marker genes
# Select the resolution for clusters
cluster_key        = "leiden_r1.5"
cluster_bname      = "leiden_r1.5"

# Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby=cluster_key, key_added='rank_genes_{0}'.format(cluster_key), n_genes=adata.shape[1])

# 4.1.2) Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_{0}'.format(cluster_key), fontsize=12, show=False)
plt.savefig("{0}/06_{1}_{2}_marker_genes_ranking.png".format(plotsDir, bname, cluster_key) , bbox_inches='tight', dpi=175); plt.close('all')

# SCSA 
# Convert the normal marker file to scsa compatible marker file 
marker_file        = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/tregCNS_markers_various_celltypes_V1_TK.txt'
convert_markerlist_to_custom_scsa_format(marker_file, species="mouse")

# 4.2.1) Annotate using SCSA a cell type annotation for single-cell RNA-seq data
marker_file        = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/tregCNS_markers_various_celltypes_V1_TK_scsa.txt'
marker_list_name   = "manec_markers_progenitor_cells_V1_scsa"
cellTypeColumnName = "CellType"
adata = annotate_with_scsa(adata, dataDir, cluster_key, marker_file, projName, marker_list_name, cellTypeColumnName, analysis_stage_num='05', analysis_stage='scsa_clusters', foldchange=2.0, pvalue=0.05)

# 4.2.2) Plot selected features
nElementsCluster = len(adata.obs['CellType'].value_counts())
colPalette       = get_color_palette(nElementsCluster)
cluster_features = ['CellType', 'sampleID', 'condition']
plot_umap_tsne(adata, plotsDir, "{0}_{1}_CellType_UMAP_TSNE".format(bname, marker_list_name), main_title = "{0} Celltype TSNE/UMAP".format(marker_list_name), features=cluster_features, analysis_stage_num='05', analysis_stage='scsa_clusters', color_palette=colPalette)

# 4.2.3) leiden UMAPs and TSNEs
plot_selected_cluster_umap_tsne(adata, plotsDir, bname, main_title = 'Leiden clustering', features=cluster_key, additional_features=cluster_features, analysis_stage_num='05', analysis_stage='scsa_clusters', color_palette=colPalette)

# 4.2.4) Plot individual samples
plot_individual_cluster_umap(adata, plotsDir, bname, cluster_key='CellType', cluster_bname='CellType', analysis_stage_num='05', analysis_stage='scsa_clusters', final_color_palette=colPalette)

# 4.2.5) Plot bar plots of cluster and samples
plot_barplots(adata, plotsDir, bname, cluster_key='CellType', cluster_bname='CellType', analysis_stage_num='05', analysis_stage='scsa_clusters', color_palette=colPalette)

# 4.3) Save the cellType assigned adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/05_SCSA_Leiden_Clustering_{1}_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/05_SCSA_Leiden_Clustering_{1}_adata.h5ad" .format(dataDir, projName); scsaleidenadata  = sc.read_h5ad(adatafile)
# adata = scsaleidenadata.copy()
