# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

ipython3 # ipython is Python 3.6.8 (default, Jan 14 2019, 11:02:34)

import scanpy as sc
from gjainPyLib import *

# Import data
# Raw data
rawData     = sc.read_text('input/manac/T_cellranger_filtered_manac_counts_raw_genesymbols.txt'); rawData.var_names_make_unique(); 
rawData
# AnnData object with n_obs × n_vars = 6398 × 30958

# Denoised data (DCA)
dcaData     = sc.read_text('output/manac/01_preprocessed/02_denoised_T_cellranger_filtered_manac_counts_raw_genesymbols.txt'); dcaData.var_names_make_unique(); 
dcaData
# AnnData object with n_obs × n_vars = 5881 × 15165 

# scTransform normalized only raw data 
normRawData = sc.read_text('output/manac/01_preprocessed/03_T_cellranger_filtered_manac_counts_normalizedRaw_genesymbols.txt'); normRawData.var_names_make_unique(); normRawData
# AnnData object with n_obs × n_vars = 6398 × 16371

# scTransform normalized denoised data 
normDcaData = sc.read_text('output/manac/01_preprocessed/03_T_cellranger_filtered_manac_counts_normalizedDCA_genesymbols.txt'); normDcaData.var_names_make_unique(); normDcaData
# AnnData object with n_obs × n_vars = 5881 × 15165

# For cell cycle QC: https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb

# Quality control 
sc.pp.filter_genes(rawData    , min_cells=50); sc.pp.filter_cells(rawData    , min_genes=2000); rawData.var_names_make_unique() 
# AnnData object with n_obs × n_vars = 5880 × 12899
#     obs: 'n_genes' 12899
#     var: 'n_cells' 5880
sc.pp.filter_genes(dcaData    , min_cells=50); sc.pp.filter_cells(dcaData    , min_genes=2000); dcaData.var_names_make_unique() 
# AnnData object with n_obs × n_vars = 5881 × 15165 
sc.pp.filter_genes(normRawData, min_cells=50); sc.pp.filter_cells(normRawData, min_genes=2000); normRawData.var_names_make_unique()
# AnnData object with n_obs × n_vars = 5776 × 12899 
sc.pp.filter_genes(normDcaData, min_cells=50); sc.pp.filter_cells(normDcaData, min_genes=2000); normDcaData.var_names_make_unique() 
# AnnData object with n_obs × n_vars = 5866 × 15165 


mito_genes = rawData.var_names.str.startswith('MT-'); rawData.obs['percent_mito'] = np.sum(rawData[:, mito_genes].X, axis=1) / np.sum(rawData.X, axis=1)
rawData.obs['n_counts'] = rawData.X.sum(axis=1);rawData = rawData[rawData.obs['n_genes'] < 2000, :]; rawData = rawData[rawData.obs['percent_mito'] < 0.05, :]

# Perform additional preprocessing
nn=10; npcs=20
sc.pp.log1p(rawData); sc.pp.scale(rawData, max_value=10); sc.tl.pca(rawData, svd_solver='arpack'); sc.pp.neighbors(rawData, n_neighbors=nn, n_pcs=npcs); sc.tl.umap(rawData)
sc.pp.log1p(dcaData); sc.pp.scale(dcaData, max_value=10); sc.tl.pca(dcaData, svd_solver='arpack'); sc.pp.neighbors(dcaData, n_neighbors=nn, n_pcs=npcs); sc.tl.umap(dcaData)
sc.pp.log1p(normRawData); sc.pp.scale(normRawData, max_value=10); sc.tl.pca(normRawData, svd_solver='arpack'); sc.pp.neighbors(normRawData, n_neighbors=nn, n_pcs=npcs); sc.tl.umap(normRawData)
sc.pp.log1p(normDcaData); sc.pp.scale(normDcaData, max_value=10); sc.tl.pca(normDcaData, svd_solver='arpack'); sc.pp.neighbors(normDcaData, n_neighbors=nn, n_pcs=npcs); sc.tl.umap(normDcaData)

output_file = 'output/manac/01_preprocessed/plots/03_cellranger_filtered_manac_data'; create_dir("{0}".format(get_file_info(output_file)[0]))
plot_data()

# Plot 
def plot_data(aedata, output_file, tag):
    ''' Plot various plots '''
    sc.pl.highest_expr_genes(aedata, n_top=10,show=False)
    plt.savefig("{0}_{1}_top_10_genes.png".format(get_file_info(output_file)[3], tag) , bbox_inches='tight'); plt.close('all')
    sc.pl.violin(aedata, ['n_genes', 'n_counts'], jitter=0.4, multi_panel=True,show=False)
    plt.savefig("{0}_{1}_violinplot_qc.png".format(get_file_info(output_file)[3], tag) , bbox_inches='tight'); plt.close('all')
    sc.pl.scatter(aedata, x='n_counts', y='n_genes',show=False)
    plt.savefig("{0}_{1}_scatterplot_ngenes_ncounts.png".format(get_file_info(output_file)[3], tag) , bbox_inches='tight'); plt.close('all')

    sc.pl.umap(aedata, color='RdBu', use_raw=False,show=False)
    plt.savefig("{0}_{1}_umap_GL1.png".format(get_file_info(output_file)[3], tag) , bbox_inches='tight'); plt.close('all')

# clit + D + D
#