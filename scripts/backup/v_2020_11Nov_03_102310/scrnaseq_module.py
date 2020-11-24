def myfunc():
    print('hello')

def  perform_qc(adata, plotsDir, bname, batch_key='sampleID', num_neighbors=15, perplexity=30, random_state=2105):
  """
  Perform QC analysis and generate QC plots

  Returns:
      [annData]: A quality controlled filtered annotated data matrix with raw counts
  """

  # Make a copy of qcadata
  qcadata = adata.copy()
  
  # 1.2.1) Calculate QC covariates
  print("- Shape {0}".format(qcadata.to_df().shape))
  sc.pp.calculate_qc_metrics(qcadata, inplace=True) # we now have many additional data types in the obs slot:
  qcadata.obs['n_counts']   = qcadata.X.sum(1)
  qcadata.obs['log_counts'] = np.log(qcadata.obs['n_counts'])
  qcadata.obs['n_genes']    = (qcadata.X > 0).sum(1)
  qcadata
  # obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'n_counts', 'log_counts', 'n_genes'
  # var: 'gene_ids', 'feature_types', 'genome', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'

  # 1.2.2) Calculate mitochondrial/ribosomal genes fraction (percentage)
  # For each cell compute fraction of counts in mito/ribo genes vs. all genes
  mt_gene_mask         = [gene.startswith('mt-') for gene in qcadata.var_names]
  qcadata.obs['mt_frac'] = qcadata.X[:, mt_gene_mask].sum(1)/qcadata.obs['n_counts']
  rb_gene_mask         = [gene.startswith(("Rps","Rpl")) for gene in qcadata.var_names]
  qcadata.obs['rb_frac'] = qcadata.X[:, rb_gene_mask].sum(1)/qcadata.obs['n_counts']

  # Plot unfiltered QC data
  print("- Plot unfiltered QC data")
  qc_plots(qcadata, plotsDir, "{0}_unfiltered".format(bname))

  # Write the information in the log file and save it in a dictionary
  f      = open("{0}/data/01_QC_info_cells_genes.txt".format(get_file_info(plotsDir)[0]), "w")
  qcdict = OrderedDict()

  # 1.2.3) Filter cells according to identified QC thresholds:
  f.write("- Unfiltered rawqcadata shape: {0}".format(qcadata.shape))
  qcdict['Unfiltered raw qc adata'] = "{0} cells and {1} genes".format(qcadata.shape[0], qcadata.shape[1])
  print('Total number of cells: {:d}'.format(qcadata.n_obs))
  f.write('\n\n- Total number of cells: {:d}'.format(qcadata.n_obs))
  qcdict['Total number of cells'] = '{:d}'.format(qcadata.n_obs)

  sc.pp.filter_cells(qcadata, min_counts = minCountPerCell)
  print('Number of cells after min count filter ({: ^6d}): {:d}'.format(minCountPerCell, qcadata.n_obs))
  f.write('\n\t- Number of cells after min count filter ({: ^6d}): {:d}'.format(minCountPerCell, qcadata.n_obs))
  qcdict['Number of cells after min count filter ({:d})'.format(minCountPerCell)] = '{:d}'.format(qcadata.n_obs)

  sc.pp.filter_cells(qcadata, max_counts = maxCountPerCell)
  print('Number of cells after max count filter ({: ^6d}): {:d}'.format(maxCountPerCell, qcadata.n_obs))
  f.write('\n\t- Number of cells after max count filter ({: ^6d}): {:d}'.format(maxCountPerCell, qcadata.n_obs))
  qcdict['Number of cells after max count filter ({:d})'.format(maxCountPerCell)] = '{:d}'.format(qcadata.n_obs)

  qcadata = qcadata[qcadata.obs['mt_frac'] < mtGenesFilter]
  print('Number of cells after MT filter   ({: ^6f}): {:d}'.format(mtGenesFilter, qcadata.n_obs))
  f.write('\n\t- Number of cells after MT filter   ({: ^6f}): {:d}'.format(mtGenesFilter, qcadata.n_obs))
  qcdict['Number of cells after MT filter ({:0.2f})'.format(mtGenesFilter)] = '{:d}'.format(qcadata.n_obs)

  qcadata = qcadata[qcadata.obs['rb_frac'] < rbGenesFilter]
  print('Number of cells after Ribo filter ({: ^6f}): {:d}'.format(rbGenesFilter, qcadata.n_obs))
  f.write('\n\t- Number of cells after Ribo filter ({: ^6f}): {:d}'.format(rbGenesFilter, qcadata.n_obs))
  qcdict['Number of cells after Ribo filter ({:0.2f})'.format(rbGenesFilter)] = '{:d}'.format(qcadata.n_obs)

  sc.pp.filter_cells(qcadata, min_genes = minGenesPerCell)
  print('Number of cells after gene filter ({: ^6d}): {:d}'.format(minGenesPerCell, qcadata.n_obs))
  f.write('\n\t- Number of cells after gene filter ({: ^6d}): {:d}'.format(minGenesPerCell, qcadata.n_obs))
  qcdict['Number of cells after gene filter ({:d})'.format(minGenesPerCell)] = '{:d}'.format(qcadata.n_obs)

  # 1.2.4) Filter genes according to identified QC thresholds:
  print('Total number of genes: {:d}'.format(qcadata.n_vars))
  f.write('\n\n- Total number of genes: {:d}'.format(qcadata.n_vars))
  qcdict['Total number of genes'] = '{:d}'.format(qcadata.n_vars)
  sc.pp.filter_genes(qcadata, min_cells=minCellsPergene)
  print('Number of genes after minCellsPergene filter ({: ^6d}): {:d}'.format(minCellsPergene, qcadata.n_vars))
  f.write('\n\t- Number of genes after minCellsPergene filter ({: ^6d}): {:d}'.format(minCellsPergene, qcadata.n_vars))
  qcdict['Number of genes after minCellsPergene filter ({:d})'.format(minCellsPergene)] = '{:d}'.format(qcadata.n_vars)

  print("- Filtered rawqcadata shape: {0}".format(qcadata.shape))
  f.write("\n\n- Filtered rawqcadata shape: {0}".format(qcadata.shape))
  qcdict['Filtered raw qc adata'] = "{0} cells and {1} genes".format(qcadata.shape[0], qcadata.shape[1])
  f.close()

  # Convert dictionary to dataframe
  # qcDF = pd.DataFrame(qcdict, columns=qcdict.keys(), index=[0])
  qcDF = pd.DataFrame(qcdict, index=[0]).T
  qcDF.style.export_png("{0}/01_{1}_QC_info_cells_genes.png".format(plotsDir, bname))

  # Plot filtered QC data
  print("- Plot filtered QC data")
  qc_plots(qcadata, plotsDir, "{0}_filtered".format(bname))

  # Calculations for the visualizations
  qcadata = calculate_umap_tsne(qcadata, num_neighbors=15, perplexity=30, random_state=2105)
  
  print("- Plot filtered QC data UMAPs and TSNEs")
  # Get the base features
  features = ['log_counts', 'mt_frac', 'rb_frac']
  # If batch key is provided then plot that as well
  if batch_key is not None:
    features.extend(batch_key)

  plot_umap_tsne(qcadata, plotsDir, "{0}_filtered_UMAP_TSNE".format(bname), features=features, analysis_stage_num='01', analysis_stage='raw', color_palette=sc.pl.palettes.vega_20_scanpy)

  return qcadata

def calculate_umap_tsne(qcadata, num_neighbors=15, perplexity=30, random_state=2105):
  """[summary] Calculations for the visualizations

  Returns:
      [type]: [description]
  """
  # Compute variable genes
  # We first need to define which features/genes are important in our dataset to distinguish cell types. 
  # For this purpose, we need to find genes that are highly variable across cells, 
  # which in turn will also provide a good separation of the cell clusters.
  sc.pp.highly_variable_genes(qcadata, flavor='cell_ranger')
  print('\n','Number of highly variable genes: {:d}'.format(np.sum(qcadata.var['highly_variable'])))

  # Compute PCA coordinates, loadings and variance decomposition
  sc.pp.pca(qcadata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = random_state)

  # Compute a neighborhood graph
  sc.pp.neighbors(qcadata, random_state = random_state, n_neighbors=num_neighbors)

  # Embed the neighborhood graph using UMAP
  sc.tl.umap(qcadata, random_state = random_state, n_components=3)

  # Embed the neighborhood graph using TSNE
  sc.tl.tsne(qcadata , random_state = random_state, n_pcs=50, perplexity=perplexity)

  return qcadata

def qc_plots(qcadata, plotsDir, bname):
  """
  Plot QC matrices

  Returns: None
  """
  n = 5
  fig = plt.figure(figsize=(40,15))
  # Sample quality plots
  ax = fig.add_subplot(2, n, 1); t1  = sc.pl.violin(qcadata, ['n_genes_by_counts', 'n_counts'], jitter=0.4, size=2, log=True, cut=0, ax = ax, show=False)
  ax = fig.add_subplot(2, n, 6); t2  = sc.pl.violin(qcadata, ['mt_frac','rb_frac'], jitter=0.4, size=2, log=False, cut=0, ax = ax, show=False)
  # 1.2.4) Thresholdingecision based on counts
  ax = fig.add_subplot(2, n, 2); p3  = sns.distplot(qcadata.obs['n_counts'], kde=False, ax = ax, bins=50); #plt.show()
  ax = fig.add_subplot(2, n, 7); p4  = sns.distplot(qcadata.obs['n_counts'][qcadata.obs['n_counts']<2000], kde=False, ax = ax, bins=50); #plt.show()
  # 1.2.5) Thresholding decision based on genes
  ax = fig.add_subplot(2, n, 3); p6  = sns.distplot(qcadata.obs['n_genes'], kde=False, ax = ax, bins=50); # plt.show()
  ax = fig.add_subplot(2, n, 8); p7  = sns.distplot(qcadata.obs['n_genes'][qcadata.obs['n_genes']<1000], kde=False, ax = ax, bins=50); # plt.show()
  # 1.2.6) mt fraction plots
  ax = fig.add_subplot(2, n, 4); p8  = sc.pl.scatter(qcadata, 'n_counts', 'n_genes', color='mt_frac', ax = ax, show=False)
  ax = fig.add_subplot(2, n, 9); p9  = sc.pl.scatter(qcadata[qcadata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', ax = ax, show=False)
  # 1.2.7) rb fraction plots
  ax = fig.add_subplot(2, n, 5); p10 = sc.pl.scatter(qcadata, 'n_counts', 'n_genes', color='rb_frac', ax = ax, show=False)
  ax = fig.add_subplot(2, n, 10);p11 = sc.pl.scatter(qcadata[qcadata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='rb_frac', ax = ax, show=False)
  # Save plot
  plt.tight_layout()
  plt.savefig("{0}/01_raw_{1}_QC_matrices.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

def plot_umap_tsne(qcadata, plotsDir, bname, main_title = 'Filtered_raw', features=None, analysis_stage_num='01', analysis_stage='raw', color_palette=sc.pl.palettes.vega_20_scanpy):
  """[summary]

  Returns:
      [type]: [description]
  """
  # 1.2.11) Plot visualizations
  sc.pl.pca_scatter(qcadata, color='n_counts',show=False)
  plt.savefig("{0}/01_raw_{1}_clustering_ncounts_PCA.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

  # Get plot parameters
  subplot_title_fontsize = 12
  subplot_title_width    = 50
  nrows                  = len(features)
  fig                    = plt.figure(figsize=(25,6*nrows))
  fig.suptitle(main_title)
  
  # Plot individual marker genes
  m=n=o=0
  for i, mfeature in enumerate(features):
    m=1+i*3; n=2+i*3; o=3+i*3;
    # print("- {0}) {4}: m={1}, n={2}, o={3}".format(i, m, n, o, mfeature))
    ax = fig.add_subplot(nrows, 3, m);                  sc.pl.tsne(qcadata, legend_loc=None     , ax=ax, color=mfeature, palette=color_palette, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}".format(mfeature),subplot_title_width)),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows, 3, n);                  sc.pl.umap(qcadata, legend_loc=None     , ax=ax, color=mfeature, palette=color_palette, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}".format(mfeature),subplot_title_width)),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows, 3, o, projection='3d'); sc.pl.umap(qcadata                      , ax=ax, color=mfeature, palette=color_palette, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("\n".join(wrap("{0}".format(mfeature),subplot_title_width)),fontsize= subplot_title_fontsize)
    
  plt.tight_layout()
  plt.savefig("{0}/{3}_{1}_{2}.png".format(plotsDir, bname, analysis_stage, analysis_stage_num) , bbox_inches='tight', dpi=175); plt.close('all')

def plot_individual_cluster_umap(qcadata, plotsDir, bname, cluster_key='sampleID', cluster_bname='sampleID', analysis_stage_num='01', analysis_stage='raw', final_color_palette=sc.pl.palettes.vega_20_scanpy):
  """[summary]

  Returns:
      [type]: [description]
  """
  # Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
  cluster_key_groups = qcadata.obs[cluster_key].cat.categories.tolist()
  cluster_cell_count = qcadata.obs[cluster_key].value_counts().to_dict()

  # # Get the color palette as variable from the string
  # # https://stackoverflow.com/questions/1373164/how-do-i-create-a-variable-number-of-variables
  # final_color_palette = getattr(sc.pl.palettes, color_palette)

  # UMAPs
  subplot_title_fontsize = 12
  subplot_title_width    = 50
  nrows  = len(cluster_key_groups) + 1
  fig = plt.figure(figsize=(30, 7*nrows))
  fig.suptitle("{0} TSNE/UMAP".format(cluster_key))
  # Main Cluster
  ax = fig.add_subplot(nrows,3, 1); sc.pl.tsne(qcadata, legend_loc=None, ax=ax, color=cluster_key, palette=final_color_palette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
  ax = fig.add_subplot(nrows,3, 2); sc.pl.umap(qcadata, legend_loc=None, ax=ax, color=cluster_key, palette=final_color_palette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
  ax = fig.add_subplot(nrows,3, 3, projection='3d'); sc.pl.umap(qcadata, ax=ax, color=cluster_key, palette=final_color_palette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
  # Partial visualizaton of a subset of groups in embedding
  # m=3; n=4
  for i,b in enumerate(cluster_key_groups):
    print(i, b)
    m=4+i*3; n=5+i*3; o=6+i*3;
    ax = fig.add_subplot(nrows,3, m);                  sc.pl.tsne(qcadata[qcadata.obs[cluster_key]== b], legend_loc=None, ax=ax, color=cluster_key, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows,3, n);                  sc.pl.umap(qcadata[qcadata.obs[cluster_key]== b], legend_loc=None, ax=ax, color=cluster_key, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows,3, o, projection='3d'); sc.pl.umap(qcadata[qcadata.obs[cluster_key]== b]                 , ax=ax, color=cluster_key, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
    m+=1; n+=1

  plt.tight_layout()
  plt.savefig("{0}/{4}_{3}_{1}_{2}_TSNE_UMAP_individual_clusters.png".format(plotsDir, bname, cluster_bname, analysis_stage, analysis_stage_num) , bbox_inches='tight', dpi=175); plt.close('all')

def save_adata_to_excel(qcadata, dataDir, outputFileName, obs_additional_colnames=None, append_new_colnames=False, obs_colnames=None, subset_genes=None):
  """
  Name:
    save_adata_to_excel
  
  Description: 
    Save the annData as a tab separated file with genes and observation columns for each cell

  Parameters:
    qcadata             (annData) : The ann data object that is used at the current stage of the analysis
    dataDir             (str)     : Path the data directory
    outputFileName      (str)     : Output filename
    obs_colnames        (list)    : Default column names that will be added with the genes
                                    Default = ['n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'log_counts', 'n_genes', 'mt_frac', 'rb_frac']
 
    obs_additional_colnames (list): New or additional column names that will be added with the genes
                                    Example: ['batch', 'sampleID', 'AllTumor_hgMycIresCd2', 'humanMyc', 'humanMycMappedToMouseMyc', 'gap', 'ires', 'humanCd2', 'louvain_r1']

    subset_genes        (list)    : A small subset of genes that should be saved.
                                    Example: ['Hbb-bs','Hba-a1','Hba-a2']
    append_new_colnames (bool)    : False
    
  Return: 
    None
  """

  # Set default column names
  # Note: It is not advisable to use mutatble default arguments to a function. 
  #       Hence the default argument is None and then assign the data type into the function. 
  #       Help: https://nikos7am.com/posts/mutable-default-arguments/
  if (obs_additional_colnames is None): obs_additional_colnames = []
  if (subset_genes is None): subset_genes = []
  if (obs_colnames is None): obs_colnames = ['n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'log_counts', 'n_genes', 'mt_frac', 'rb_frac']

  # Extend the column names to the default column names
  if (append_new_colnames):
    obs_colnames.extend(obs_additional_colnames)

  # print(obs_colnames)

  # Get the dataframes for count data and observations
  countsDF = qcadata.to_df()
  obsDF    = qcadata.obs.copy()

  # Subset the observation dataframe if the parameters are set
  subset_obsDF = obsDF[obs_colnames].copy()

  # Subset the counts dataframe if the parameters are set
  if(subset_genes):
    qcadata_expressed_gene_list = qcadata.var.index.tolist()
    subset_expressed_gene_list  = [x for x in subset_genes if x in qcadata_expressed_gene_list]
    subset_countsDF             = countsDF[subset_expressed_gene_list].copy()
  else:
    subset_countsDF             = countsDF.copy()

  # Merge the selected observation and count dataframes
  print("\n- Merge the selected observation and count dataframes")
  merged_subset_obs_countDF = pd.concat([subset_obsDF, subset_countsDF], axis=1)

  # Get the descriptive statistics for the DF
  merged_subset_obs_countDF = merged_subset_obs_countDF.describe()[1:].append(merged_subset_obs_countDF)
  
  # Save the dataframe
  merged_subset_obs_countDF.T.to_csv("{0}/{1}".format(dataDir, outputFileName), sep='\t', header=True, index=True, index_label="cellIDs")
  
########################
# Normalization Modules
########################
def cpm_norm(adata, plotsDir, bname):
  """ 
  # Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells
  """
  cpmadata = adata.copy()
  sc.pp.normalize_total(cpmadata, target_sum=1e4)
  # Logarithmize the data
  sc.pp.log1p(cpmadata) 
  cpmadata.raw = cpmadata
  return cpmadata

def scran_norm(adata, plotsDir, bname):
  """
  Normalization using SCRAN

  Args:
      adata ([anndata]): [description]
  """
  scranadata = adata.copy()
  # Perform a clustering for scran normalization in clusters
  adata_pp = scranadata.copy()
  sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
  sc.pp.log1p(adata_pp)
  sc.pp.pca(adata_pp, n_comps=15)
  sc.pp.neighbors(adata_pp)
  sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)

  # Preprocess variables for scran normalization
  input_groups = adata_pp.obs['groups']
  data_mat = scranadata.X.T
  
  # Run scran in R
  # https://rpy.sourceforge.io/rpy2/doc-dev/html/introduction.html
  robjects.globalenv["data_mat"] = data_mat
  robjects.globalenv["input_groups"] = input_groups
  size_factors = scran.computeSumFactors(data_mat, clusters=input_groups)
  
  # Delete adata_pp
  del adata_pp
  
  # Visualize the estimated size factors
  scranadata.obs['size_factors'] = size_factors
  fig = plt.figure(figsize=(16,6))
  fig.suptitle('Estimated size factors')
  ax = fig.add_subplot(1, 2, 1)
  sc.pl.scatter(scranadata, 'size_factors', 'n_counts', ax=ax, show=False)
  ax = fig.add_subplot(1, 2, 2)
  sc.pl.scatter(scranadata, 'size_factors', 'n_genes', ax=ax, show=False)
  plt.tight_layout()
  plt.savefig("{0}/02_norm_{1}_scran_sizefactors_plots.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  # Keep the count data in a counts layer
  scranadata.layers["counts"] = scranadata.X.copy()
  # Normalize adata 
  scranadata.X /= scranadata.obs['size_factors'].values[:,None]
  sc.pp.log1p(scranadata)
  # Store the full data set in 'raw' as log-normalised data for statistical testing
  scranadata.raw = scranadata

  return scranadata

########################
# Batch Correction Modules
########################
def combat_bc(adata, plotsDir, bname, batchkey='tissueID'):
  """
  Batch correction using combat

  Args:
      adata ([anndata]): [description]
  """
  combatadata = adata.copy()
  sc.pp.combat(combatadata, key=batchkey)
  return combatadata

def scanorama_bc(qcadata, plotsDir, bname, batchkey='batch'):
  """
  Batch correction using scanorama

  Args:
      adata ([anndata]): [description]
  """
  qcadata2 = sc.AnnData(X=qcadata.X, var=qcadata.var, obs = qcadata.obs)
  # Variable genes for the full dataset
  sc.pp.highly_variable_genes(qcadata2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = batchkey)
  var_genes_batch = qcadata2.var.highly_variable_nbatches > 0
  var_select = qcadata2.var.highly_variable_nbatches > 1
  var_genes = var_select.index[var_select]
  # Split per batch into new objects.
  # batches = ['0','1','2','3','4','5','6','7','8','9','10','11','12']
  batches = adata.obs[batchkey].cat.categories.tolist()
  qcalladata = {}
  for batch in batches:
      qcalladata[batch] = qcadata2[qcadata2.obs[batchkey] == batch,]

  # Subset the individual dataset to the same variable genes as in MNN-correct.
  qcalladata2 = dict()
  for ds in qcalladata.keys():
      print(ds)
      qcalladata2[ds] = qcalladata[ds][:,var_genes]
  # Convert to list of AnnData objects
  qcnormadatas = list(qcalladata2.values())
  # Run scanorama.integrate
  qcscanorama  = scanorama.integrate_scanpy(qcnormadatas, dimred = 50,)
  # Make into one matrix.
  qcall_s = np.concatenate(qcscanorama)
  print(qcall_s.shape)

  return qcall_s

def cell_cycle_correction(adata, plotsDir, bname):
  """

  """
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
  plt.savefig("{0}/02_norm_{1}_scran_cell_cycle_plots.png".format(plotsDir, bname) , bbox_inches='tight', dpi=150); plt.close('all')


###############################
# MARKER GENES MODULES
###############################

def import_marker_genes_list(marker_file, species="mouse"):
  """
  Name:
    import_marker_genes_list
  
  Description: 
    Read the marker genes into a pandas dataframe and return a marker genes dictionary

  Parameters:
    marker_file (Str):  - Marker genes file name with path.
                        - Expected format: CellTypes <tab> MarkerGenes (comma separated list) with the header CellTypes and MarkerGenes

                        ┌────────────────────────────┬────────────────────────────────────────────────────────────────────────┐
                        │         CellTypes          │                              MarkerGenes                               │
                        ├────────────────────────────┼────────────────────────────────────────────────────────────────────────┤
                        │ Pit_cell                   │ Gm26917,Gm42418,Gm23935                                                │
                        │ Epithelial_cell_MT4_high   │ Krt10,Krt1,Mt4,Krtdap,Dapl1,Fabp5,Dmkn,Lypd3,Lgals7,Plin2              │
                        │ Epithelial_cell_Clca1_high │ Krt20,Clca1,Psapl1,Lypd8,Fcgbp,Phgr1,Lgals4,Fabp2,Sptssb,2210407C18Rik │
                        └────────────────────────────┴────────────────────────────────────────────────────────────────────────┘
                       
                       - Note: Please avoid using spaces or other non-ascii characters. Underscores are okay. 

  Return: 
    marker_genes_dict (Dict): A dictionary containing marker genes list. The keys are celltypes and values are marker genes
  """
  marker_genesDF    = pd.read_csv(marker_file, sep="\t", header=0, index_col=None)
  if species == 'mouse':
    marker_genes_dict = marker_genesDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.lower().capitalize() for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()
  elif species == 'human':
    marker_genes_dict = marker_genesDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.upper()              for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()
  else:
    marker_genes_dict = marker_genesDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x for x in n.split(',')] for i in g.values.tolist() for n in i]))).to_dict()
  return marker_genes_dict

def plot_manual_marker_list_genes(adata, markerDir, bname, cluster_key, genespresent, marker_genes_cellTypes, marker_list_name):
  """
  Generate the UMAPs and TSNEs for each marker categories

  Args:
      adata ([anndata]): [description]

  Desc:
      # Get genes that are present in the adata
      # l = ['a', 'b', 'c', 'd', 'f']
      # d = {'A':['a','x','c'], 'B':['c','d'],'C':['x','y']}

      # for k,v in d.items():
      #   nv = [x for x in v if x in l]
      #   d[k] = nv

  """
  for k,v in marker_genes_cellTypes.items():
    print("\n- Original list {0}: {1}".format(k,v))
    validgenes = [x for x in v if x in genespresent]
    ids = np.in1d(adata.var_names,validgenes)
    print("- Genes present {0}: {1}".format(k,validgenes))

    subplot_title_fontsize = 12
    subplot_title_width    = 50
    ngenes = len(validgenes)
    nrows  = ngenes + 2
    adata.obs['{0}_marker_expr'.format(k)] = adata.X[:,ids].mean(1)

    fig = plt.figure(figsize=(25,6*nrows))
    fig.suptitle(marker_list_name)
    # Plot cluster
    ax = fig.add_subplot(nrows, 3, 1                 ); sc.pl.tsne(adata, legend_loc='on data', ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3,                  show=False, title="{0} TSNE".format(cluster_key))
    ax = fig.add_subplot(nrows, 3, 2);                  sc.pl.umap(adata, legend_loc='on data', ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format(cluster_key))
    ax = fig.add_subplot(nrows, 3, 3, projection='3d'); sc.pl.umap(adata                      , ax=ax, color="{0}".format(cluster_key), palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format(cluster_key))
    
    # Plots mean marker genes
    ax = fig.add_subplot(nrows, 3, 4);                  sc.pl.tsne(adata, legend_loc=None     , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))
    ax = fig.add_subplot(nrows, 3, 5);                  sc.pl.umap(adata, legend_loc=None     , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))
    ax = fig.add_subplot(nrows, 3, 6, projection='3d'); sc.pl.umap(adata                      , ax=ax, color='{0}_marker_expr'.format(k), color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("Mean {0}".format("\n".join(wrap("{0}:{1}".format(k,validgenes),subplot_title_width)),fontsize= subplot_title_fontsize))
    
    # Plot individual marker genes
    m=n=o=0
    for i, mgene in enumerate(validgenes):
      m=7+i*3; n=8+i*3; o=9+i*3;
      # print("- {0}) {4}: m={1}, n={2}, o={3}".format(i, m, n, o, mgene))
      ax = fig.add_subplot(nrows, 3, m);                  sc.pl.tsne(adata, legend_loc=None     , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
      ax = fig.add_subplot(nrows, 3, n);                  sc.pl.umap(adata, legend_loc=None     , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
      ax = fig.add_subplot(nrows, 3, o, projection='3d'); sc.pl.umap(adata                      , ax=ax, color=mgene, color_map=mymap, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("\n".join(wrap("{0}:{1}".format(k,mgene),subplot_title_width)),fontsize= subplot_title_fontsize)
      
    plt.tight_layout()
    plt.savefig("{0}/{1}_{2}_{3}_TSNE_UMAP.png".format(markerDir, bname, marker_list_name, k) , bbox_inches='tight', dpi=100); plt.close('all')

def plot_additional_marker_gene_visualization(adata, markerDir, bname, cluster_key, genespresent, marker_genes_dict, marker_list_name, analysis_stage_num='04', analysis_stage='norm'):
  """
  Generate the UMAPs and TSNEs for each marker categories

  Args:
      adata ([anndata]): [description]

  Desc:
      # Get genes that are present in the adata
      # l = ['a', 'b', 'c', 'd', 'f']
      # d = {'A':['a','x','c'], 'B':['c','d'],'C':['x','y']}

      # for k,v in d.items():
      #   nv = [x for x in v if x in l]
      #   d[k] = nv

  """
  # 5.2) Other marker gene visualization
  # 5.2.1) Get clean marker dict
  adata_expressed_genes = adata.var.index.tolist()
  marker_genes_filtered_dict = defaultdict()
  for k,v in marker_genes_dict.items():
    new_genes_list = [x for x in v if x in adata_expressed_genes]
    if new_genes_list:
      marker_genes_filtered_dict[k] = new_genes_list

  # Calculate the dendrogram again
  sc.tl.dendrogram(adata, groupby=cluster_key)

  # 5.2.2) Dot plots: 
  #   - The dotplot visualization provides a compact way of showing per group, 
  #     the fraction of cells expressing a gene (dot size) and the mean expression 
  #     of the gene in those cell (color scale)
  #   - The use of the dotplot is only meaningful when the counts matrix contains 
  #     zeros representing no gene counts. dotplot visualization does not work for 
  #     scaled or corrected matrices in which cero counts had been replaced by other values.
  sc.pl.dotplot(adata, marker_genes_filtered_dict, groupby=cluster_key, log=True, figsize=(40,12), show=False, dendrogram=True)
  plt.savefig("{0}/{1}_{2}_{3}_{4}_{5}_dotplot.png".format(plotsDir, analysis_stage_num, analysis_stage, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')

  # 5.2.3) Matrix plots: 
  #   - The matrixplot shows the mean expression of a gene in a group by 
  #     category as a heatmap. 
  #   - In contrast to dotplot, the matrix plot can be used with corrected 
  #     and/or scaled counts. By default raw counts are used.
  # sc.pl.matrixplot(adata, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False,cmap='Reds',  figsize=(40,12), standard_scale='group', show=False)
  sc.pl.matrixplot(adata, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False,cmap='Reds',  figsize=(40,12), standard_scale='group', show=False)
  plt.savefig("{0}/{1}_{2}_{3}_{4}_{5}_scaled_matrixplot.png".format(plotsDir, analysis_stage_num, analysis_stage, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')
  sc.pl.matrixplot(adata, marker_genes_filtered_dict, groupby=cluster_key, dendrogram=True, use_raw=False, cmap='Reds', figsize=(40,12), standard_scale='group', vmin=0.5, show=False)
  plt.savefig("{0}/{1}_{2}_{3}_{4}_{5}_scaled_vmin0_05_matrixplot.png".format(plotsDir, analysis_stage_num, analysis_stage, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')

  # 5.2.4) Tracksplots: 
  #   - The track plot shows the same information as the heatmap, but, 
  #     instead of a color scale, the gene expression is represented by height
  sc.pl.tracksplot(adata, marker_genes_filtered_dict, groupby=cluster_key, log=True, dendrogram=True, show=False, figsize=(50,30))
  plt.savefig("{0}/{1}_{2}_{3}_{4}_{5}_trackplot.png".format(plotsDir, analysis_stage_num, analysis_stage, bname, cluster_key, marker_list_name) , bbox_inches='tight', dpi=175); plt.close('all')


def plot_barplots(adata, plotsDir, bname, cluster_key='sampleID', cluster_bname='sampleID', analysis_stage_num='01', analysis_stage='raw', color_palette="vega_20"):
  """
  Plot separate bar plots, coloured in by cluster annotation, for each tissue
  """
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
  plt.savefig("{0}/{4}_{3}_{1}_{2}_tissueID_cluster_barplot.png".format(plotsDir, bname, cluster_bname, analysis_stage, analysis_stage_num) , bbox_inches='tight', dpi=175); plt.close('all')

# 