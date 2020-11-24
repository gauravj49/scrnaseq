# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)

import scanpy as sc
# mp.use('Agg') # to use matplotlib without X11
from gjainPyLib import *


projName        = "stomach1001_mouse"
fileType        = ""
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/02_exploratory_analysis/{0}".format(projName); create_dir("{0}".format(output_dir))
minGenesPerCell = 200
minCellsPergene = 3

# Import data
# inputFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/01_raw_T_{0}_cellranger_filtered_manec_counts_genesymbols.txt'.format(projName)
# inputFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/02_denoised_DCA_T_{0}_cellranger_filtered_manec_counts_genesymbols.txt'.format(projName)
# inputFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/03_normalizedRawCounts_T_{0}_cellranger_filtered_manec_counts_genesymbols.txt'.format(projName)
inputFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/03_normalizedDCACounts_T_{0}_cellranger_filtered_manec_counts_genesymbols.txt'.format(projName)

bname = get_file_info(inputFile)[1]

# Get the AnnData object
aedata = sc.read_text(inputFile); aedata.var_names_make_unique(); aedata
# AnnData object with n_obs (cells) × n_vars (genes) = 4254 × 13502 

# Calculate qc metrics for visualization
# aedata.obs.columns                                                                                                                                                           Index(['n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts',
    #    'log1p_total_counts', 'pct_counts_in_top_50_genes',
    #    'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes',
    #    'pct_counts_in_top_500_genes'],
    #   dtype='object')


# 2 Pre-processing and visualization

# # 2.1 Quality control
# # Data quality control can be split into cell QC and gene QC. Typical quality measures for assessing the quality of a cell include the number of molecule counts (UMIs), the number of expressed genes, and the fraction of counts that are mitochondrial. A high fraction of mitochondrial reads being picked up can indicate cell stress, as there is a low proportion of nuclear mRNA in the cell. It should be noted that high mitochondrial RNA fractions can also be biological signals indicating elevated respiration

# obsPlotFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/02_exploratory_analysis/01_raw_T_bulk997_mouse_cellranger_filtered_manec_counts_genesymbols_Obs_QC_metrices.png'; create_dir(get_file_info(obsPlotFile)[0])

# sc.pp.calculate_qc_metrics(aedata, inplace=True)
# # sns.jointplot("log1p_total_counts", "log1p_n_genes_by_counts", data=aedata.obs, kind='regplot')
# sns.scatterplot(x="n_genes_by_counts", y="total_counts", data=aedata.obs)
# # sc.pl.violin(aedata.obs, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True)
# plt.savefig(obsPlotFile)

# # Quality control - calculate QC covariates
# aedata.obs['n_counts'] = aedata.X.sum(1)
# aedata.obs['log_counts'] = np.log(aedata.obs['n_counts'])
# aedata.obs['n_genes'] = (aedata.X > 0).sum(1)

# mt_gene_mask = [gene.startswith('mt-') for gene in aedata.var_names]
# aedata.obs['mt_frac'] = aedata.X[:, mt_gene_mask].sum(1)/aedata.obs['n_counts']

# # Quality control - plot QC metrics
# # Sample quality plots
# t1 = sc.pl.violin(aedata, 'n_counts', size=2, log=True, cut=0)
# t2 = sc.pl.violin(aedata, 'mt_frac')


# # Data quality summary plots
# p1 = sc.pl.scatter(aedata, 'n_counts', 'n_genes', color='mt_frac')
# p2 = sc.pl.scatter(aedata[aedata.obs['n_counts']<10000], 'n_counts', 'n_genes', color='mt_frac')

# For cell cycle QC: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
# For batch effects: https://nbviewer.jupyter.org/github/Teichlab/bbknn/blob/master/examples/pancreas.ipynb

# Basic filtering
sc.pp.filter_cells(adata, min_genes=minGenesPerCell)
sc.pp.filter_genes(adata, min_cells=minGenesPerCell)
adata.var_names_make_unique()


# Plot top 10 genes
sc.pl.highest_expr_genes(adata, n_top=10, show = False)
plt.savefig("{0}/{1}_top_10_genes.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')

# Mitochondrial genes
# for each cell compute fraction of counts in mito genes vs. all genes
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1)

# A violin plot of the computed quality measures.
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, show=False)
plt.savefig("{0}/{1}_violinplot_qc.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')

# Scatter plot of ngenes_ncounts
# Ideal should be on the 45degree line
sc.pl.scatter(adata, x='n_counts', y='n_genes', show = False)
plt.savefig("{0}/{1}_scatterplot_ncounts_ngenes.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.scatter(adata, x='n_counts', y='percent_mito', show=False)
plt.savefig("{0}/{1}_scatterplot_ncounts_pcMito.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')

# Actual filtering and TPM normalize
adata = adata[adata.obs['n_genes'] < 1000, :]
adata = adata[adata.obs['percent_mito'] < 0.25, :]
# Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells.
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

# Logarithmize the data.
sc.pp.log1p(adata)

# Set the .raw attribute of AnnData object to the logarithmized raw gene expression for later use in differential testing and visualizations of gene expression
adata.raw = adata

# Identify highly-variable genes.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=1, min_disp=0.5)
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig("{0}/{1}_highly_variable_genes.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')



