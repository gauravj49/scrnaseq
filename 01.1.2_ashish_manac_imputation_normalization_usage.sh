# pwd

# ############### Imputation ########################

ipython3

# Denoising using DCA
import scanpy as sc
from dca.api import dca
from gjainPyLib import *

ipython # Python 3.6.8 (default, Jan 14 2019, 11:02:34)
aedata = sc.read_text('input/manac/T_cellranger_filtered_manac_counts_raw_genesymbols.txt')
aedata.var_names_make_unique()
aedata
# Out[30]: AnnData object with n_obs × n_vars = 6398 × 30958


output_file = 'output/manac/01_imputation/denoised'
create_dir("{0}".format(get_file_info(output_file)[0]))
sc.pp.filter_genes(aedata, min_cells=10)
sc.pp.filter_cells(aedata, min_genes=2000)
aedata.var_names_make_unique()
aedata
# AnnData object with n_obs × n_vars = 5881 × 15165 
#     obs: 'n_genes' 15165
#     var: 'n_cells' 5881

# Denoise the data on raw counts
dca(aedata)
# DCA: Successfully preprocessed 15165 genes and 5881 cells.

# Save the denoised data
aedata.to_df().to_csv("{0}_T_cellranger_filtered_manac_counts_raw_genesymbols.txt".format(get_file_info(output_file)[3]), sep='\t', header=True, index=True, index_label="CellId")

aedata.to_df().T.to_csv("{0}_cellranger_filtered_manac_counts_raw_genesymbols.txt".format(get_file_info(output_file)[3]), sep='\t', header=True, index=True, index_label="GeneSymbol")

# clit + D + D


 
# ############### Normalization ########################

