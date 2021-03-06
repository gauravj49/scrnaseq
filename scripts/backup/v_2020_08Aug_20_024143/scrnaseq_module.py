def myfunc():
    print('hello')

def  perform_qc(adata, plotsDir, bname):
  """
  Perform QC analysis and generate QC plots
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

  # Write the information in the log file
  f = open("{0}/data/01_QC_info_cells_genes.txt".format(get_file_info(plotsDir)[0]), "w")
  # 1.2.3) Filter cells according to identified QC thresholds:
  f.write("- Unfiltered rawqcadata shape: {0}".format(qcadata.shape))
  print('Total number of cells: {:d}'.format(qcadata.n_obs))
  f.write('\n\n- Total number of cells: {:d}'.format(qcadata.n_obs))
  sc.pp.filter_cells(qcadata, min_counts = minCountPerCell)
  print('Number of cells after min count filter: {:d}'.format(qcadata.n_obs))
  f.write('\n\t- Number of cells after min count filter: {:d}'.format(qcadata.n_obs))

  sc.pp.filter_cells(qcadata, max_counts = maxCountPerCell)
  print('Number of cells after max count filter: {:d}'.format(qcadata.n_obs))
  f.write('\n\t- Number of cells after max count filter: {:d}'.format(qcadata.n_obs))

  qcadata = qcadata[qcadata.obs['mt_frac'] < mtGenesFilter]
  print('Number of cells after MT filter  : {:d}'.format(qcadata.n_obs))
  f.write('\n\t- Number of cells after MT filter  : {:d}'.format(qcadata.n_obs))
  qcadata = qcadata[qcadata.obs['rb_frac'] < rbGenesFilter]
  print('Number of cells after Ribo filter: {:d}'.format(qcadata.n_obs))
  f.write('\n\t- Number of cells after Ribo filter: {:d}'.format(qcadata.n_obs))

  sc.pp.filter_cells(qcadata, min_genes = minGenesPerCell)
  print('Number of cells after gene filter: {:d}'.format(qcadata.n_obs))
  f.write('\n\t- Number of cells after gene filter: {:d}'.format(qcadata.n_obs))

  # 1.2.4) Filter genes according to identified QC thresholds:
  print('Total number of genes: {:d}'.format(qcadata.n_vars))
  f.write('\n\n- Total number of genes: {:d}'.format(qcadata.n_vars))
  sc.pp.filter_genes(qcadata, min_cells=minCellsPergene)
  print('Number of genes after minCellsPergene filter: {:d}'.format(qcadata.n_vars))
  f.write('\n\t- Number of genes after minCellsPergene filter: {:d}'.format(qcadata.n_vars))

  print("- Filtered rawqcadata shape: {0}".format(qcadata.shape))
  f.write("\n\n- Filtered rawqcadata shape: {0}".format(qcadata.shape))
  f.close()

  # Plot filtered QC data
  print("- Plot filtered QC data")
  qc_plots(qcadata, plotsDir, "{0}_filtered".format(bname))
  
  print("- Plot filtered QC data UMAPs")
  plot_raw_umap(qcadata, plotsDir, "{0}_filtered".format(bname), 25)

  return qcadata

def qc_plots(qcadata, plotsDir, bname):
  """
  Plot QC metrics
  """
  n = 5
  fig = plt.figure(figsize=(15,40))
  # Sample quality plots
  ax = fig.add_subplot(n, 2, 1); t1  = sc.pl.violin(qcadata, ['n_genes_by_counts', 'n_counts'], jitter=0.4, size=2, log=True, cut=0, ax = ax, show=False)
  ax = fig.add_subplot(n, 2, 2); t2  = sc.pl.violin(qcadata, ['mt_frac','rb_frac'], jitter=0.4, size=2, log=False, cut=0, ax = ax, show=False)
  # 1.2.4) Thresholdingecision based on counts
  ax = fig.add_subplot(n, 2, 3); p3  = sns.distplot(qcadata.obs['n_counts'], kde=False, ax = ax, bins=50); #plt.show()
  ax = fig.add_subplot(n, 2, 4); p4  = sns.distplot(qcadata.obs['n_counts'][qcadata.obs['n_counts']<2000], kde=False, ax = ax, bins=50); #plt.show()
  # 1.2.5) Thresholding decision based on genes
  ax = fig.add_subplot(n, 2, 5); p6  = sns.distplot(qcadata.obs['n_genes'], kde=False, ax = ax, bins=50); # plt.show()
  ax = fig.add_subplot(n, 2, 6); p7  = sns.distplot(qcadata.obs['n_genes'][qcadata.obs['n_genes']<1000], kde=False, ax = ax, bins=50); # plt.show()
  # 1.2.6) mt fraction plots
  ax = fig.add_subplot(n, 2, 7); p8  = sc.pl.scatter(qcadata, 'n_counts', 'n_genes', color='mt_frac', ax = ax, show=False)
  ax = fig.add_subplot(n, 2, 8); p9  = sc.pl.scatter(qcadata[qcadata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='mt_frac', ax = ax, show=False)
  # 1.2.7) rb fraction plots
  ax = fig.add_subplot(n, 2, 9); p10 = sc.pl.scatter(qcadata, 'n_counts', 'n_genes', color='rb_frac', ax = ax, show=False)
  ax = fig.add_subplot(n, 2, 10);p11 = sc.pl.scatter(qcadata[qcadata.obs['n_counts']<2000], 'n_counts', 'n_genes', color='rb_frac', ax = ax, show=False)
  # Save plot
  plt.tight_layout()
  plt.savefig("{0}/01_raw_{1}_QC_matrices.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

def plot_raw_umap(qcadata, plotsDir, bname, num_neighbors=15):
  # # For debugging purpose
  # qcadata=rawadata.copy()
  # 1.2.9) Compute variable genes
  # We first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.
  sc.pp.highly_variable_genes(qcadata, flavor='cell_ranger')
  print('\n','Number of highly variable genes: {:d}'.format(np.sum(qcadata.var['highly_variable'])))

  # 1.2.10) Calculations for the visualizations
  sc.pp.pca(qcadata, n_comps=50, use_highly_variable=True, svd_solver='arpack', random_state = 2105)
  sc.pp.neighbors(qcadata, random_state = 2105, n_neighbors=num_neighbors)
  sc.tl.umap(qcadata, random_state = 2105, n_components=3)

  # 1.2.11) Plot visualizations
  sc.pl.pca_scatter(qcadata, color='n_counts',show=False)
  plt.savefig("{0}/01_raw_{1}_clustering_ncounts_PCA.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')
  # UMAPS
  fig = plt.figure(figsize=(20,32))
  # 2D projection
  ax = fig.add_subplot(4, 2, 1);                  sc.pl.umap(qcadata   ,                  ax=ax, color='sampleID'  , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format('log_counts'))
  ax = fig.add_subplot(4, 2, 3);                  sc.pl.umap(qcadata   ,                  ax=ax, color='log_counts', palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="{0} UMAP".format('log_counts'))
  ax = fig.add_subplot(4, 2, 5);                  sc.pl.umap(qcadata   , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="mt_frac UMAP")
  ax = fig.add_subplot(4, 2, 7);                  sc.pl.umap(qcadata   , legend_loc=None, ax=ax, color="rb_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False, title="rb_frac UMAP")
  # 3D projection
  ax = fig.add_subplot(4, 2, 2, projection='3d'); sc.pl.umap(qcadata   , legend_loc=None, ax=ax, color='sampleID'  , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format('log_counts'))
  ax = fig.add_subplot(4, 2, 4, projection='3d'); sc.pl.umap(qcadata   , legend_loc=None, ax=ax, color='log_counts', palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="{0} UMAP".format('log_counts'))
  ax = fig.add_subplot(4, 2, 6, projection='3d'); sc.pl.umap(qcadata   , legend_loc=None, ax=ax, color="mt_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="mt_frac UMAP")
  ax = fig.add_subplot(4, 2, 8, projection='3d'); sc.pl.umap(qcadata   , legend_loc=None, ax=ax, color="rb_frac"   , palette=sc.pl.palettes.vega_20, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False, title="rb_frac UMAP")
  plt.tight_layout()
  plt.savefig("{0}/01_raw_{1}_sampleID_logCounts_mt_rb_frac_UMAP.png".format(plotsDir, bname) , bbox_inches='tight', dpi=175); plt.close('all')

def plot_individual_cluster_umap(qcadata, plotsDir, bname, cluster_key='sampleID', cluster_bname='sampleID', analysis_stage_num='01', analysis_stage='raw', color_palette="vega_20"):
  """[summary]

  Returns:
      [type]: [description]
  """
  # Get number of groups for the cluster_key (cluster_key_groups,number_of_cells)
  cluster_key_groups = qcadata.obs[cluster_key].cat.categories.tolist()
  cluster_cell_count = qcadata.obs[cluster_key].value_counts().to_dict()

  # Get the color palette as variable from the string
  # https://stackoverflow.com/questions/1373164/how-do-i-create-a-variable-number-of-variables
  final_color_palette = getattr(sc.pl.palettes, color_palette)
  # Louvain UMAPs
  subplot_title_fontsize = 12
  subplot_title_width    = 50
  ncols  = len(cluster_key_groups) + 1
  fig = plt.figure(figsize=(20, 7*ncols))
  fig.suptitle("{0} UMAP".format(cluster_key))
  # Main Louvain Cluster
  ax = fig.add_subplot(ncols,2, 1); sc.pl.umap(qcadata, legend_loc=None, ax=ax, color=cluster_key, palette=final_color_palette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False)
  ax = fig.add_subplot(ncols,2, 2, projection='3d'); sc.pl.umap(qcadata, ax=ax, color=cluster_key, palette=final_color_palette, size=100, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False)
  # Partial visualizaton of a subset of groups in embedding
  m=3; n=4
  for i,b in enumerate(cluster_key_groups):
    print(i, b)
    # ax = fig.add_subplot(ncols,2, i+m);                  sc.pl.umap(qcadata, legend_loc=None, ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
    # ax = fig.add_subplot(ncols,2, i+n, projection='3d'); sc.pl.umap(qcadata                 , ax=ax, color=cluster_key, groups=[b], size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(ncols,2, i+m);                  sc.pl.umap(qcadata[qcadata.obs[cluster_key]== b], legend_loc=None, ax=ax, color=cluster_key, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, show=False);  ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
    ax = fig.add_subplot(ncols,2, i+n, projection='3d'); sc.pl.umap(qcadata[qcadata.obs[cluster_key]== b]                 , ax=ax, color=cluster_key, size=50, edgecolor='k', linewidth=0.05, alpha=0.9, hspace=0.35, wspace=0.3, projection='3d', show=False); ax.set_title("{0}: {1} cells".format(b, cluster_cell_count[b]),fontsize= subplot_title_fontsize)
    m+=1; n+=1

  plt.tight_layout()
  plt.savefig("{0}/{4}_{3}_{1}_{2}_UMAP_individual_clusters.png".format(plotsDir, bname, cluster_bname, analysis_stage, analysis_stage_num) , bbox_inches='tight', dpi=175); plt.close('all')

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
  batches = ['0','1','2','3','4','5','6','7','8','9','10','11','12']
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
