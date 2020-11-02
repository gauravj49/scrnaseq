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
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/{0}/allSamples_ohneS515".format(projName); create_dir("{0}".format(output_dir))
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
adatas          = [sc.read_10x_mtx(f) for f in tissueFilenames]
adatas

# 1.2) Get the dictionary of tissue ids
# tissueIdDict =  {'0' :'S503', '1' :'S504', '2' :'S505', '3' :'S508', '4' :'S509', '5' :'S511', '6' :'S512', '7' :'S514', '8' :'S515', '9' :'S516', '10':'S517', '11':'S518', '12':'S519'}
tissueIdDict =  {'0' :'S503', '1' :'S504', '2' :'S505', '3' :'S508', '4' :'S509', '5' :'S511', '6' :'S512', '7' :'S514', '8' :'S516', '9' :'S517', '10':'S518', '11':'S519'}

# 1.3) Merge 10x datasets for different mices
adata = adatas[0].concatenate(adatas[1:])

# 1.4) Make variable names unique
adata.var_names_make_unique()

# Convert the sparse count matrices to dense represntation
adata.X = adata.X.toarray()

# 1.5) Add tissue id column for the batches
adata.obs['tissueID'] = adata.obs['batch'].map(tissueIdDict)
# In [16]: adata.obs
# Out[16]:
#                     batch tissueID
# AAACAAACAGCTATGA-0      0     S503
# AAACAAACCTACGAGC-0      0     S503

# 1.6) Calculate and plot QC covariates/metrices
rawadata = perform_qc(adata, plotsDir, bname)
# - Unfiltered rawqcadata shape: (28632, 55471)
#     Total number of cells: 28632
#     Number of cells after min count filter: 26139
#     Number of cells after max count filter: 26116
#     Number of cells after MT filter  : 26098
#     Number of cells after Ribo filter: 25329
#     Number of cells after gene filter: 24891

#     Total number of genes: 55471
#     Number of genes after minCellsPergene filter: 16617
# - Filtered rawqcadata shape: (24891, 16617)

# 1.7) Save the filtered raw adata into a file
# Write the adata object to file
adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata.write(adatafile)
# # Read back the filtered raw adata object
# adatafile  = "{0}/01_raw_{1}_adata.h5ad" .format(dataDir, projName); rawadata  = sc.read_h5ad(adatafile)

########################
# 2) Normalization
########################
cpmadata   = cpm_norm(rawadata, plotsDir, bname)
scranadata = scran_norm(rawadata, plotsDir, bname)

# ########################
# # 3) Cell cycle correction
# ########################
# cell_cycle_correction()

# 3.2) Save the normalized adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/02_norm_cpm_adata.h5ad" .format(dataDir, projName); cpmadata.write(adatafile)
adatafile  = "{0}/02_norm_scran_adata.h5ad" .format(dataDir, projName); scranadata.write(adatafile)
# # Read back the corrected adata object
# adatafile  = "{0}/02_norm_cpm_adata.h5ad"  .format(dataDir, projName); cpmadata   = sc.read_h5ad(adatafile)
# adatafile  = "{0}/02_norm_scran_adata.h5ad".format(dataDir, projName); scranadata = sc.read_h5ad(adatafile)
# normadata = adata.copy()

########################
# 4) Batch correction
########################
# cpmcombatadata = combat_bc(cpmadata, plotsDir, bname, batchkey='tissueID')
# scranscanoramaadata = scanorama_bc(scranadata, plotsDir, bname, batchkey='tissueID'):

normadata = scranadata.copy()
rawadatas = adatas.copy()
adata     = normadata.copy()

# 5) Technical correction: Batch Correction using Scanorama
# 5.1) Detect variable genes
# As the stored AnnData object contains scaled data based on variable genes, we need to make a new object with the raw counts and normalized it again. Variable gene selection should not be performed on the scaled data object, only do normalization and log transformation before variable genes selection.
scranadata2 = sc.AnnData(X=scranadata.X, var=scranadata.var, obs = scranadata.obs)
#variable genes for the full dataset
sc.pp.highly_variable_genes(scranadata2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'batch')
var_genes_batch = scranadata2.var.highly_variable_nbatches > 0
var_select = scranadata2.var.highly_variable_nbatches > 1
var_genes = var_select.index[var_select]
# Split per batch into new objects.
# batches = ['0','1','2','3','4','5','6','7','8','9','10','11','12']
batches = list(tissueIdDict.keys())
scranalldata = {}
for batch in batches:
    scranalldata[batch] = scranadata2[scranadata2.obs['batch'] == batch,]
# Subset the individual dataset to the same variable genes as in MNN-correct.
scranalldata2 = dict()
for ds in scranalldata.keys():
    print(ds)
    scranalldata2[ds] = scranalldata[ds][:,var_genes]
# Convert to list of AnnData objects
comnormadatas = list(scranalldata2.values())
# Run scanorama.integrate
scranscanorama  = scanorama.integrate_scanpy(comnormadatas, dimred = 50,)
# Make into one matrix.
scranall_s = np.concatenate(scranscanorama)
print(scranall_s.shape)
# Add to the AnnData object
scranscanoramaadata = adata.copy()
scranscanoramaadata.obsm["SC"] = scranall_s
# Calculations for the visualizations
sc.pp.highly_variable_genes(scranscanoramaadata, flavor='cell_ranger', n_top_genes=4000)
sc.pp.pca(scranscanoramaadata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
sc.pp.neighbors(scranscanoramaadata, random_state = 2105, use_rep = "SC")
sc.tl.umap(scranscanoramaadata, random_state = 2105, n_components=3)

# Calculations for the visualizations
sc.pp.neighbors(rawadata, random_state = 2105)
sc.tl.umap(rawadata, random_state = 2105, n_components=3)

adata = scranscanoramaadata.copy()

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

# 3.2) Save the normalized batch corrected adata into a file
# Write the adata and cadata object to file
adatafile  = "{0}/03_norm_all_batchCorrection_adata.h5ad" .format(dataDir, projName); adata.write(adatafile)
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
