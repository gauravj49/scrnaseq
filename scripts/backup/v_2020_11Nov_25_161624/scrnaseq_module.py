def myfunc():
    print('hello')

def import_datasets(input_h5_files, sampleIdDict=None, batch_key=None):
  """
  Name:
    import_datasets
  
  Description: 
    Read 10x H5 datasets and return the annotation data object

  Parameters:
    input_h5_files (str)    : Input H5 file or list of files. For example:
                              - input_h5_files  = 'input/manec/stomachBaseline/3520_Antrum_filtered_feature_bc_matrix.h5'
                              - or
                              - input_h5_files  = ['input/manec/stomachBaseline/3306_tumorA_filtered_feature_bc_matrix.h5',
                                                   'input/manec/stomachBaseline/3306_tumorB_filtered_feature_bc_matrix.h5'
                                                  ]
    sampleIdDict   (dict)   : Dictionary of sample/tissue ids. This should be in exact same order as the input files. 
                              For example:
                              sampleIdDict =  {
                                              '0':'3306_tumorA',
                                              '1':'3306_tumorB',
                                              }
    batch_key      (str)    : Name of the batch key that will be used during the downstream analysis. For example: sampleID

  Return: 
    adata          (annData): The ann data object that is used at the current stage of the analysis
  """

  if isinstance(input_h5_files, str):
    # 1.1) Reading the data in the anndata object individually
    adata   = sc.read_10x_h5(input_h5_files)
  elif isinstance(input_h5_files, list):
    # 1.1.1) Reading the data in the anndata object individually
    adatas  = [sc.read_10x_h5(f) for f in input_h5_files]
    
    # 1.1.2) Make variable names unique
    for ad in adatas: print(ad.shape); ad.var_names_make_unique()

    # 1.1.4) Merge 10x datasets for different mices
    adata = adatas[0].concatenate(adatas[1:])

    # 1.1.5) Add tissue id column for the batches
    adata.obs[batch_key]  = adata.obs['batch'].map(sampleIdDict)
    # In [16]: adata.obs
    # Out[16]:
    #                     batch sampleID
    # AAACAAACAGCTATGA-0      0     S503
    # AAACAAACCTACGAGC-0      0     S503

  # 1.2) Make variable names unique
  adata.var_names_make_unique()

  # 1.3) Convert the sparse count matrices to dense representation
  adata.X = adata.X.toarray()

  # Return the ann data object
  return adata


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
    features.append(batch_key)
    print(features)

  plot_umap_tsne(qcadata, plotsDir, "{0}_filtered_UMAP_TSNE".format(bname), features=features, analysis_stage_num='01', analysis_stage='raw', color_palette=sc.pl.palettes.vega_20_scanpy)

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

def calculate_umap_tsne(qcadata, num_neighbors=15, perplexity=30, early_exaggeration=12, random_state=2105):
  """[summary] Calculations for the visualizations

  perplexity : float, int (default: 30)
               The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. 
               Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. 
               The choice is not extremely critical since t-SNE is quite insensitive to this parameter.

  early_exaggeration : Union[float, int] (default: 12)
              Controls how tight natural clusters in the original space are in the embedded space and how much space will be between them. 
              For larger values, the space between natural clusters will be larger in the embedded space. Again, the choice of this 
              parameter is not very critical. If the cost function increases during initial optimization, the early exaggeration factor or 
              the learning rate might be too high.
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
  sc.tl.tsne(qcadata , random_state = random_state, n_pcs=50, perplexity=perplexity, early_exaggeration=early_exaggeration)

  return qcadata

def plot_umap_tsne(qcadata, plotsDir, bname, main_title = 'Filtered_raw', features=None, analysis_stage_num='01', analysis_stage='raw', color_palette=sc.pl.palettes.vega_20_scanpy):
  """[summary]

  Returns:
      [type]: [description]
  """
  # 1.2.11) Plot visualizations
  sc.pl.pca_scatter(qcadata, color='n_counts',show=False)
  plt.savefig("{0}/{3}_{1}_{2}_ncounts_PCA.png".format(plotsDir, bname, analysis_stage, analysis_stage_num), bbox_inches='tight', dpi=175); plt.close('all')

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
  for i,b in enumerate(cluster_key_groups):
    print(i, b)
    m=4+i*3; n=5+i*3; o=6+i*3;
    ax = fig.add_subplot(nrows,3, m);                  sc.pl.tsne(qcadata[qcadata.obs[cluster_key]== b], legend_loc=None, ax=ax, color=cluster_key, size=300, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows,3, n);                  sc.pl.umap(qcadata[qcadata.obs[cluster_key]== b], legend_loc=None, ax=ax, color=cluster_key, size=300, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows,3, o, projection='3d'); sc.pl.umap(qcadata[qcadata.obs[cluster_key]== b]                 , ax=ax, color=cluster_key, size=300, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
    m+=1; n+=1

  plt.tight_layout()
  plt.savefig("{0}/{4}_{1}_{2}_{3}_TSNE_UMAP_individual_clusters.png".format(plotsDir, bname, cluster_bname, analysis_stage, analysis_stage_num) , bbox_inches='tight', dpi=175); plt.close('all')

def plot_selected_cluster_umap_tsne(qcadata, plotsDir, bname, main_title = 'Filtered_raw', features='leiden_r1', additional_features=None, analysis_stage_num='01', analysis_stage='raw', color_palette=sc.pl.palettes.vega_20_scanpy):
  """[summary]

  feature             = cluster_key
  additional_features = ['sampleID', 'condition']
  Returns:
      [type]: [description]
  """
  # Convert the feature string to list
  if isinstance(features, str):
    features = features.split(" ")

  # Add additional features to the features list
  if (additional_features is None): additional_features = []
  if additional_features:
    if isinstance(additional_features, str):
      # Append the additional feature to the features list
      features.append(additional_features)
    elif isinstance(additional_features, list):
      # Extend the column names to the default features list
      features.extend(additional_features)

  # Get plot parameters
  subplot_title_fontsize = 12
  subplot_title_width    = 50
  nrows                  = len(features)
  fig                    = plt.figure(figsize=(50,8*nrows))
  fig.suptitle(main_title)
  
  # Plot leiden/louvain UMAPs and TSNEs
  m=n=o=p=q=0
  for i, mfeature in enumerate(features):
    m=1+i*5; n=2+i*5; o=3+i*5; p=4+i*5; q=5+i*5 
    # print("- {0}) {6}: m={1}, n={2}, o={3}, p={4}, q={5}".format(i, m, n, o, p, q, mfeature))
    ax = fig.add_subplot(nrows, 5, m);                  sc.pl.tsne(qcadata, legend_loc='on data', ax=ax, color=mfeature, palette=color_palette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}".format(mfeature),subplot_title_width)),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows, 5, n);                  sc.pl.tsne(qcadata, legend_loc=None     , ax=ax, color=mfeature, palette=color_palette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}".format(mfeature),subplot_title_width)),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows, 5, o);                  sc.pl.umap(qcadata, legend_loc='on data', ax=ax, color=mfeature, palette=color_palette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}".format(mfeature),subplot_title_width)),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows, 5, p);                  sc.pl.umap(qcadata, legend_loc=None     , ax=ax, color=mfeature, palette=color_palette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("\n".join(wrap("{0}".format(mfeature),subplot_title_width)),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(nrows, 5, q, projection='3d'); sc.pl.umap(qcadata                      , ax=ax, color=mfeature, palette=color_palette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("\n".join(wrap("{0}".format(mfeature),subplot_title_width)),fontsize= subplot_title_fontsize)
    
  # Save the plot
  plt.tight_layout()
  plt.savefig("{0}/{3}_{1}_{2}_2D3D_UMAP_TSNE.png".format(plotsDir, bname, analysis_stage, analysis_stage_num) , bbox_inches='tight', dpi=175); plt.close('all')


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

########################
# CLUSTERING MODULES
########################
def get_color_palette(nelements=15):
  """
  Name:
    get_color_palette
  
  Description: 
    - Get a color palette based on number of clusters or samples
    - Provides a list of colors to use for plotting categorical annotation groups

  Parameters:
    nelements  (str) : Number of elements (ex. number of clusters or samples)

  Return: 
    colPalette (list): The color palette based on number of clusters or samples
  """
  # Default colors
  

  if nelements <= 20:
    colPalette = sc.pl.palettes.vega_20_scanpy
  elif 21 <= nelements <= 28:
    # https://epub.wu.ac.at/1692/1/document.pdf
    colPalette = sc.pl.palettes.zeileis_28
  else:
    # http://godsnotwheregodsnot.blogspot.de/2012/09/color-distribution-methodology.html
    colPalette = sc.pl.palettes.godsnot_102
  return colPalette

def calculate_plot_clustering(qcadata, plotsDir, bname, main_title = 'Clustering all resolution TSNE/UMAP', additional_features=None, analysis_stage_num='03', analysis_stage='clustering', color_palette=sc.pl.palettes.vega_20_scanpy, clustering_algorithm='leiden', random_state = 2105):
  """
  Name:
    calculate_plot_clustering
  
  Description: 
    Calculate and visualize the clustering using highly variable genes

  Parameters:
    qcadata             (annData) : The ann data object that is used at the current stage of the analysis
    plotsDir            (str)     : Path the plots directory
    bname               (str)     : Output file basename
    
    
  Return: 
    qcadata             (annData) : The ann data object that is used at the current stage of the analysis
  """
  # Check if correct clustering algorithm is passed. If not, then default to leiden
  if clustering_algorithm not in ['leiden', 'lovain']:
    print("WARNING: Could not found clustering algorithm.\n\t- Values allowed: leiden or louvain.\n\t- Defaulting to leiden.")
    clustering_algorithm = 'leiden'

  # Get resolution list
  cluster_resolution = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2]
  # Perform clustering - using highly variable genes
  print("- Perform clustering - using highly variable genes")
  for cres in cluster_resolution:
    res_key = '{0}_r{1}'.format(clustering_algorithm, cres)
    if clustering_algorithm == 'leiden':
      sc.tl.leiden(adata, resolution=cres, key_added=res_key, random_state=random_state)
    elif clustering_algorithm == 'louvain':
      sc.tl.louvain(adata, resolution=cres, key_added=res_key, random_state=random_state)
    print(adata.obs[res_key].value_counts())

  # Visualize the clustering and how this is reflected by different technical covariates
  print("- Visualize the clustering and how this is reflected by different technical covariates")
  cluster_features = ['{0}_r{1}'.format(clustering_algorithm, cres) for cres in cluster_resolution]
  plot_umap_tsne(adata, plotsDir, "{0}_clustering_all_leiden_UMAP_TSNE".format(bname), main_title = 'Clustering all resolution TSNE/UMAP', features=cluster_features, analysis_stage_num=analysis_stage_num, analysis_stage=analysis_stage, color_palette=color_palette)

  return qcadata

def annotate_with_scsa(qcadata, dataDir, cluster_key, marker_file, projName, marker_list_name, cellTypeColumnName="cellType", analysis_stage_num='05', analysis_stage='scsa_clusters', foldchange=2.0, pvalue=0.05):
  """
  Name:
    annotate_with_scsa
  
  Description: 
    Run SCSA algorith on custom marker genes file and then assign annoations of 
    the clusters obtained form scsa algorithm back to the original ann data
    object. Bascially, add a new column of cellTypes that is assigned to each cluster. 

  Parameters:
    qcadata             (annData) : The ann data object that is used at the 
                                    current stage of the analysis
    
    dataDir             (str)     : Path the data directory

    projName            (str)     : Name of the project
    
    marker_list_name    (str)     : Name of the custom marker list
    
    cluster_key         (str)     : Name of the cluster key that will be used 
                                    during the downstream analysis. 
                                    Example: 'leiden_r1'

    marker_file         (Str)     : - Marker genes file name with path.
                                    - Output format: CellTypes <tab> MarkerGenes (every gene in new line) and no headers
                                      ┌─────────────────────┬─────────┐
                                      │ CD34+               │ THY1    │
                                      │ CD34+               │ ENG     │
                                      │ CD34+               │ KIT     │
                                      │ CD34+               │ PROM1   │
                                      │ Natural killer cell │ NCAM1   │
                                      │ Natural killer cell │ FCGR3A  │
                                      │ Monocytes           │ CD14    │
                                      │ Monocytes           │ FCGR1A  │
                                      │ Monocytes           │ CD68    │
                                      │ Monocytes           │ S100A12 │
                                      └─────────────────────┴─────────┘

    cellTypeColumnName  (str)     : Name of the column that will be added to the
                                    anndata object
                                    Default = 'cellType'

    analysis_stage_num  (int)     : Stage of the analysis that will be used in the
                                    output filenames. 
                                    Default: 05

    analysis_stage      (str)     : Name of the stage of the analysis that will 
                                    be used in the output filenames. 
                                    Default: 'scsa_clusters'

    foldchange          (float)   : Fold change threshold for marker filtering. 
                                    Default = 2.0

    pvalue              (float)   : P-value threshold for marker filtering.
                                    Default = 0.05
                      
  Return: 
    qcadata             (annData) : The ann data object that is used at the current stage of the analysis
  """
  # Get the scsa input data in the dataframe object
  result = qcadata.uns['rank_genes_{0}'.format(cluster_key)]
  groups = result['names'].dtype.names
  scsaDF = pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'logfoldchanges','scores','pvals']})

  # Save the scsaDF to the input scsa file 
  scsa_input_marker_ann_file = "{0}/04.2_{1}_{2}_input.csv" .format(dataDir, projName, marker_list_name)
  scsaDF.to_csv(scsa_input_marker_ann_file)

  # 4.6.3) Run scsa and get the output in the output file
  scsa_output_file = "{0}/04.3_{1}_{2}_output.csv" .format(dataDir, projName, marker_list_name)
  # os.system("python3 scripts/SCSA/SCSA.py -d scripts/SCSA/whole.db -i {0} -g Mouse -s scanpy -E --Gensymbol -f1.5 -p 0.01 -o {1} -m txt -M {2} -N --norefdb".format(scsa_input_marker_ann_file, scsa_output_file,marker_file))
  print("python3 scripts/SCSA/SCSA.py -d scripts/SCSA/whole.db -i {0} -g Mouse -s scanpy -E --Gensymbol -f{3} -p {4} -o {1} -m txt -M {2} ".format(scsa_input_marker_ann_file, scsa_output_file, marker_file, foldchange, pvalue))
  os.system("python3 scripts/SCSA/SCSA.py -d scripts/SCSA/whole.db -i {0} -g Mouse -s scanpy -E --Gensymbol -f{3} -p {4} -o {1} -m txt -M {2} ".format(scsa_input_marker_ann_file, scsa_output_file,marker_file, foldchange, pvalue))
  # Generte the excel file and do not print anything while generating the file
  os.system("python3 scripts/SCSA/SCSA.py -d scripts/SCSA/whole.db -i {0} -g Mouse -s scanpy -E --Gensymbol -f{3} -p {4} -o {1} -m ms-excel -M {2} -b --noprint".format(scsa_input_marker_ann_file, scsa_output_file,marker_file, foldchange, pvalue))

  # 4.6.4) Read into the scsa annotation file
  scsaAnnotationFile = "{0}_SCSA_Annotation.txt".format(get_file_info(scsa_output_file)[3])
  scsaadata = assign_scsa_annotaion_to_anndata(qcadata, cluster_key, scsaAnnotationFile, cellTypeColumnName="CellType", analysis_stage_num='05', analysis_stage='scsa_clusters')

  return scsaadata

def assign_scsa_annotaion_to_anndata(qcadata, cluster_key, scsaAnnotationFile, cellTypeColumnName="cellType", analysis_stage_num='05', analysis_stage='scsa_clusters'):
  """
  Name:
    assign_scsa_annotaion_to_anndata
  
  Description: 
    Assign annoations of the clusters obtained form scsa algorithm back to the 
    original ann data object. Bascially, add a new column of cellTypes that 
    is assigned to each cluster. 

  Parameters:
    qcadata             (annData) : The ann data object that is used at the 
                                    current stage of the analysis

    cluster_key         (str)     : Name of the cluster key that will be used 
                                    during the downstream analysis. 
                                    Example: 'leiden_r1'

    scsaAnnotationFile  (str)     : The output file after SCSA annotation. Example:
                                    +---------+------+--------------+-----------+-------+
                                    | Cluster | Type |   Celltype   |   Score   | Times |
                                    +---------+------+--------------+-----------+-------+
                                    |       0 | ?    | Tcell        | 1.62|1.62 |   1.0 |
                                    |       1 | ?    | Tcell|Immune | 1.62|1.62 |   1.0 |
                                    |       2 | ?    | Pit|Tuft     | 1.26|0.81 |  1.55 |
                                    |       3 | N    | -            | -         |     - |
                                    |       4 | Good | G-cells      | 2.20      |  2.79 |
                                    +---------+------+--------------+-----------+-------+

    cellTypeColumnName  (str)     : Name of the column that will be added to the
                                    anndata object

    analysis_stage_num  (int)     : Stage of the analysis that will be used in the
                                    output filenames. Example: 03

    analysis_stage      (str)     : Name of the stage of the analysis that will 
                                    be used in the output filenames. 
                                    Example: 'clustering'
                      
  Return: 
    qcadata             (annData) : The ann data object that is used at the current stage of the analysis
  """

  # Read into the scsa annotation file
  # +---------+------+-------------------------------------------+---------------------------------------+--------------------+
  # | Cluster | Type |                 Celltype                  |                 Score                 |       Times        |
  # +---------+------+-------------------------------------------+---------------------------------------+--------------------+
  # |       0 | ?    | Progenitor_at_Neck|Progenitor_cell        | 1.620185174601965|1.620185174601965   |                1.0 |
  # |       1 | ?    | Progenitor_at_Neck|Progenitor_cell        | 1.620185174601965|1.620185174601965   |                1.0 |
  # |       2 | Good | Parietal_Progenitor                       | 2.2088380249674486                    |   2.79861265843033 |
  # |       3 | ?    | Progenitor_at_Neck|Progenitor_cell        | 1.620185174601965|1.620185174601965   |                1.0 |
  # |       4 | ?    | Progenitor_at_Neck|Progenitor_cell        | 1.6201851746019649|1.6201851746019649 |                1.0 |
  # |       5 | ?    | Progenitor_Procr_high|Parietal_Progenitor | 1.1031533135536862|0.8761049595722625 | 1.2591565673732457 |
  # |       6 | ?    | Progenitor_cell|Progenitor_at_Neck        | 1.8267896507049208|1.3394332546624488 | 1.3638526924324355 |
  # |       7 | N    | -                                         | -                                     |                  - |
  # +---------+------+-------------------------------------------+---------------------------------------+--------------------+
  # Get unique elements of a list mylist = [*{*mylist}] 

  scsaAnnotationDF   = pd.read_csv(scsaAnnotationFile, sep="\t", header=0, index_col=None)
  scsaAnnotationDF.sort_values(['Cluster'], inplace=True)

  # Get a new cell type column from the annotation of the cluster_key
  qcadata.obs[cellTypeColumnName] = adata.obs[cluster_key]

  # Get the dictionary of ClusterID and CellTypes
  clusterDict = dict(zip(scsaAnnotationDF['Cluster'], scsaAnnotationDF['Celltype']))
  # Replace '-' with the clusterX 
  clusterDict = {key:"Cluster{0}".format(key) if value == '-' else value for key, value in clusterDict.items()}

  # Get a new subcluster column. There are three steps
  # 1) Convert the cateogry column to int column
  qcadata.obs[cellTypeColumnName] = qcadata.obs[cellTypeColumnName].astype(int)
  # 2) Map the dictionary to the clusters
  qcadata.obs[cellTypeColumnName]=qcadata.obs[cellTypeColumnName].map(clusterDict)
  # 3) Convert the column back to category column
  qcadata.obs[cellTypeColumnName] = qcadata.obs[cellTypeColumnName].astype('category')

  # List new categories
  qcadata.obs[cellTypeColumnName].cat.categories

  return qcadata

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
    marker_genes_dict = marker_genesDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.replace(" ", "").lower().capitalize() for x in n.rstrip(',').split(',')] for i in g.values.tolist() for n in i]))).to_dict()
  elif species == 'human':
    marker_genes_dict = marker_genesDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.replace(" ", "").upper()              for x in n.rstrip(',').split(',')] for i in g.values.tolist() for n in i]))).to_dict()
  else:
    marker_genes_dict = marker_genesDF.groupby('CellTypes')[['MarkerGenes']].apply(lambda g: list(itertools.chain.from_iterable([[x.replace(" ", "") for x in n.rstrip(',').split(',')] for i in g.values.tolist() for n in i]))).to_dict()
  return marker_genes_dict

def convert_markerlist_to_custom_scsa_format(marker_file, species="mouse"):
  """
  Name:
    convert_markerlist_to_custom_scsa_format
  
  Description: 
    Read the marker genes into a pandas dataframe and save it in a new file with SCSA custom marker gene list format

  Parameters:
    marker_file (Str):  - Marker genes file name with path.
                        - Output format: CellTypes <tab> MarkerGenes (every gene in new line) and no headers
                          ┌─────────────────────┬─────────┐
                          │ CD34+               │ THY1    │
                          │ CD34+               │ ENG     │
                          │ CD34+               │ KIT     │
                          │ CD34+               │ PROM1   │
                          │ Natural killer cell │ NCAM1   │
                          │ Natural killer cell │ FCGR3A  │
                          │ Monocytes           │ CD14    │
                          │ Monocytes           │ FCGR1A  │
                          │ Monocytes           │ CD68    │
                          │ Monocytes           │ S100A12 │
                          └─────────────────────┴─────────┘

  Return: 
    None: Saves the with marker genes file name with path
  """
  marker_genes_dict = import_marker_genes_list(marker_file, species="mouse")
  output_scsa_file  = "{0}_scsa.txt".format(get_file_info(marker_file)[3])
  with open(output_scsa_file, 'w') as f:
    for key, value in marker_genes_dict.items():
      for v in value:
        # Save the key <tab> value in the output file
        f.write('{0}\t{1}\n'.format(key, v))

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

  # Calculate the dendrogram again
  sc.tl.dendrogram(adata, groupby=cluster_key)
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
