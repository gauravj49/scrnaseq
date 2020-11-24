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

  # 1.2.3) Filter cells according to identified QC thresholds:
  print('Total number of cells: {:d}'.format(qcadata.n_obs))
  sc.pp.filter_cells(qcadata, min_counts = minCountPerCell)
  print('Number of cells after min count filter: {:d}'.format(qcadata.n_obs))

  sc.pp.filter_cells(qcadata, max_counts = maxCountPerCell)
  print('Number of cells after max count filter: {:d}'.format(qcadata.n_obs))

  qcadata = qcadata[qcadata.obs['mt_frac'] < mtGenesFilter]
  print('Number of cells after MT filter  : {:d}'.format(qcadata.n_obs))
  qcadata = qcadata[qcadata.obs['rb_frac'] < rbGenesFilter]
  print('Number of cells after Ribo filter: {:d}'.format(qcadata.n_obs))

  sc.pp.filter_cells(qcadata, min_genes = minGenesPerCell)
  print('Number of cells after gene filter: {:d}'.format(qcadata.n_obs))

  # 1.2.4) Filter genes according to identified QC thresholds:
  print('Total number of genes: {:d}'.format(qcadata.n_vars))
  sc.pp.filter_genes(qcadata, min_cells=minCellsPergene)
  print('Number of genes after minCellsPergene filter: {:d}'.format(qcadata.n_vars))

  print("- Filtered rawqcadata shape: {0}".format(qcadata.shape))

  # Plot filtered QC data
  print("- Plot filtered QC data")
  qc_plots(qcadata, plotsDir, "{0}_filtered".format(bname))

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
  size_factors = scran.computeSumFactors(data_mat, clusters=input_groups, min_mean=0.1)
  
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

def scanorama_bc(adata, plotsDir, bname, batchkey='batch'):
  """
  Batch correction using scanorama

  Args:
      adata ([anndata]): [description]
  """
  cpmadata2 = sc.AnnData(X=cpmadata.X, var=cpmadata.var, obs = cpmadata.obs)
  # Variable genes for the full dataset
  sc.pp.highly_variable_genes(cpmadata2, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = batchkey)
  var_genes_batch = cpmadata2.var.highly_variable_nbatches > 0
  var_select = cpmadata2.var.highly_variable_nbatches > 1
  var_genes = var_select.index[var_select]
  # Split per batch into new objects.
  batches = ['0','1','2','3','4','5','6','7','8','9','10','11','12']
  cpmalldata = {}
  for batch in batches:
      cpmalldata[batch] = cpmadata2[cpmadata2.obs[batchkey] == batch,]

  # Subset the individual dataset to the same variable genes as in MNN-correct.
  cpmalldata2 = dict()
  for ds in cpmalldata.keys():
      print(ds)
      cpmalldata2[ds] = cpmalldata[ds][:,var_genes]
  # Convert to list of AnnData objects
  comnormadatas = list(cpmalldata2.values())
  # Run scanorama.integrate
  cpmscanorama  = scanorama.integrate_scanpy(comnormadatas, dimred = 50,)
  # Make into one matrix.
  cpmall_s = np.concatenate(cpmscanorama)
  print(cpmall_s.shape)
  # Add to the AnnData object
  cpmscanoramaadata = adata.copy()
  cpmscanoramaadata.obsm["SC"] = cpmall_s