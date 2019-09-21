# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

ipython # ipython is Python 3.6.8 (default, Jan 14 2019, 11:02:34)

import scanpy as sc
# mp.use('Agg') # to use matplotlib without X11
from gjainPyLib import *

# Import data
# inputFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/03_normalizedDCACounts_T_bulk997_mouse_cellranger_filtered_manec_counts_genesymbols.txt'
inputFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/01_raw_T_bulk997_mouse_cellranger_filtered_manec_counts_genesymbols.txt'

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

# 2.1 Quality control
# Data quality control can be split into cell QC and gene QC. Typical quality measures for assessing the quality of a cell include the number of molecule counts (UMIs), the number of expressed genes, and the fraction of counts that are mitochondrial. A high fraction of mitochondrial reads being picked up can indicate cell stress, as there is a low proportion of nuclear mRNA in the cell. It should be noted that high mitochondrial RNA fractions can also be biological signals indicating elevated respiration

obsPlotFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/02_exploratory_analysis/01_raw_T_bulk997_mouse_cellranger_filtered_manec_counts_genesymbols_Obs_QC_metrices.png'; create_dir(get_file_info(obsPlotFile)[0])

sc.pp.calculate_qc_metrics(aedata, inplace=True)
# sns.jointplot("log1p_total_counts", "log1p_n_genes_by_counts", data=aedata.obs, kind='regplot')
sns.scatterplot(x="n_genes_by_counts", y="total_counts", data=aedata.obs)
# sc.pl.violin(aedata.obs, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True)
plt.savefig(obsPlotFile)

# Quality control - calculate QC covariates
aedata.obs['n_counts'] = aedata.X.sum(1)
aedata.obs['log_counts'] = np.log(aedata.obs['n_counts'])
aedata.obs['n_genes'] = (aedata.X > 0).sum(1)

mt_gene_mask = [gene.startswith('mt-') for gene in aedata.var_names]
aedata.obs['mt_frac'] = aedata.X[:, mt_gene_mask].sum(1)/aedata.obs['n_counts']

# Quality control - plot QC metrics
# Sample quality plots
t1 = sc.pl.violin(aedata, 'n_counts', size=2, log=True, cut=0)
t2 = sc.pl.violin(aedata, 'mt_frac')


# Data quality summary plots
p1 = sc.pl.scatter(aedata, 'n_counts', 'n_genes', color='mt_frac')
p2 = sc.pl.scatter(aedata[aedata.obs['n_counts']<10000], 'n_counts', 'n_genes', color='mt_frac')

# For cell cycle QC: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
