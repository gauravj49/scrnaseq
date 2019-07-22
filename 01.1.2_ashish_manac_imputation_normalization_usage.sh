# pwd

# ############### Imputation ########################

# Denoising using DCA
import scanpy as sc
from dca.api import dca
from gjainPyLib import *

ipython # Python 3.6.8 (default, Jan 14 2019, 11:02:34)
aedata = sc.read_text('input/manac/cellranger_filtered_manac_counts_raw_genesymbols.txt')
aedata.var_names_make_unique()
aedata
# Out[22]: AnnData object with n_obs × n_vars = 30958 × 6398


output_file = 'output/manac/inputation/01_denoised'
create_dir("{0}".format(get_file_info(output_file)[0]))
sc.pp.filter_cells(aedata, min_genes=200)
sc.pp.filter_genes(aedata, min_cells=2)
aedata.var_names_make_unique()
aedata
# AnnData object with n_obs × n_vars = 11213 × 6398 
#     obs: 'n_genes'
#     var: 'n_cells'

# Denoise the data on raw counts
dca(aedata)

# Save the denoised data
aedata.to_df().to_csv("{0}_T_cellranger_filtered_manac_counts_raw_genesymbols.txt".format(get_file_info(output_file)[3]), sep='\t', header=True, index=True, index_label="CellId")

aedata.to_df().T.to_csv("{0}_cellranger_filtered_manac_counts_raw_genesymbols.txt".format(get_file_info(output_file)[3]), sep='\t', header=True, index=True, index_label="GeneSymbol")


# ############### Normalization ########################

