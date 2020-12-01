# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

ipython # Python 3.7.0 (default, Jun 28 2018, 13:15:42)

import scanpy as sc
# mp.use('Agg') # to use matplotlib without X11
from gjainPyLib import *


projName        = "bulk997_mouse" # MANEC_merged_except1079_hMYC_forcecells
fileType        = ""
output_dir      = "/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/02_exploratory_analysis/{0}".format(projName); create_dir("{0}".format(output_dir))
minGenesPerCell = 100
minGenesPerCell = 1

# Import data
# inputFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/01_raw_T_{0}_cellranger_filtered_manec_counts_genesymbols.txt'.format(projName)
# inputFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/02_denoised_DCA_T_{0}_cellranger_filtered_manec_counts_genesymbols.txt'.format(projName)
# inputFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/03_normalizedRawCounts_T_{0}_cellranger_filtered_manec_counts_genesymbols.txt'.format(projName)
inputFile = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/output/manec/pilot2/01_preprocessed/03_normalizedDCACounts_T_{0}_cellranger_filtered_manec_counts_genesymbols.txt'.format(projName)

# Get the AnnData object
adata = sc.read_text(inputFile); adata.var_names_make_unique(); adata
# AnnData object with n_obs (cells) × n_vars (genes) = 4254 × 13502 
bname = get_file_info(inputFile)[1]

# adata = sc.read_10x_mtx(
#     '/media/rad/HDD2/temp_manec/{0}/outs/filtered_feature_bc_matrix/'.format(projName),  # the directory with the `.mtx` file
#     var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
#     cache=True)                                # write a cache file for faster subsequent reading
# adata.var_names_make_unique()
# adata
# bname = '{0}_cellranger_filtered_manec_counts_genesymbols'.format(projName)

# Convert gene names to upper
adata.var_names = adata.var_names.str.upper()

# # Split the index to get clusters
# # ┌────────────────────┬────────────┬─────────┐
# # │                    │    louvain │ mouseid │
# # ├────────────────────┼────────────┼─────────┤
# # │ AAACCCACACCGTGGT.1 │         12 │       1 │
# # │ AAACGCTAGGGTACAC.1 │         10 │       1 │
# # │ AAACGCTGTACTAGCT.1 │          8 │       1 │
# # └────────────────────┴────────────┴─────────┘
# adata.obs['mouseid'] = adata.obs.index.to_series().str.split("\.",expand=True)[1]

# # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
# sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

# # Scale each gene to unit variance. Clip values exceeding standard deviation 10.
# sc.pp.scale(adata, max_value=10)

# Principal component analysis
sc.tl.pca(adata, svd_solver='arpack')
# sc.pl.pca_variance_ratio(adata, log=True, show=False)
# plt.savefig("{0}/{1}_pca_variace_ratio.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
# results_file="{0}/{1}_pca.h5ad".format(output_dir, bname)
# adata.write(results_file)

# sc.pl.pca(adata, color='MYC', show=False)
# plt.savefig("{0}/{1}_PCA_MYC.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')


# Embedding the neighborhood graph
sc.tl.umap(adata)

# Clustering the neighborhood graph
sc.tl.louvain(adata)

# Import Markers list
marker_file = '/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/docs/stomach_marker_list_V1.txt' 
markersDF   = pd.read_csv(marker_file, sep="\t")

# loop through the markers list and generate the umaps
for nt in markersDF.itertuples(index=False):
	cellLines   = nt[0]
	cellTypes   = nt[1]
	markerGenes = nt[2].split(',')
	markerGenes.extend(['louvain'])
	#print(cellLines, cellTypes, markerGenes)

	# Plot the umaps
	try:
		sc.pl.umap(adata, color=markerGenes, use_raw=False, color_map='plasma', show=False)
		plt.savefig("{0}/{1}_UMAP_{2}_marker_genes.png".format(output_dir, bname, cellTypes) , bbox_inches='tight'); plt.close('all')
	except:
		print("\n- NOTE: {0} is not present in the dataset.".format(markerGenes))

sc.pl.umap(adata, color=['louvain','APOE', 'VWF', 'DCN', 'CD3D','CD74', 'H2-AB1', 'PF4'], use_raw=False, color_map='plasma', show=False)
plt.savefig("{0}/{1}_UMAP_top_marker_genes.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')


# sc.pl.umap(adata, color=['louvain','mouseid','GAMT','CHGA', 'CD2', 'MS4A1', 'CD3D', 'EPCAM', 'SYP', 'CD163', 'VWF', 'ATP4A', 'MUC6', 'VIM'], use_raw=False, color_map='plasma', show=False)
sc.pl.umap(adata, color=['louvain','GAMT','CD2', 'MS4A1', 'CD3D', 'EPCAM',  'VWF', 'ATP4A', 'VIM', 'CD19'], use_raw=False, color_map='plasma', show=False)
plt.savefig("{0}/{1}_UMAP_Tumorcells_marker_genes.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.umap(adata, color=['louvain','OLFM4','SOX2','LGR5','CCKBR','LRIG1','TNFRSF19'], use_raw=False, color_map='plasma', show=False)
plt.savefig("{0}/{1}_UMAP_StemCells_marker_genes.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.umap(adata, color=['OLFM4','SOX2','LGR5','CCKBR','LRIG1','TNFRSF19'], use_raw=False, color_map='plasma', show=False)
plt.savefig("{0}/{1}_UMAP_StemCells_marker_genes.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.umap(adata, color=['ID1','MBD1','RSPO2'], use_raw=False, color_map='plasma', show=False)
plt.savefig("{0}/{1}_UMAP_StomachPitCell_marker_genes.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.umap(adata, color=['LRIG1','AXIN2','CD44','ACTC1'], use_raw=False, color_map='plasma', show=False)
plt.savefig("{0}/{1}_UMAP_StomachPitProginitor_marker_genes.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')

sc.pl.umap(adata, color=['louvain','CHGA','CHGB','TAC1','TPH1','NEUROG3','SYP'], use_raw=False, color_map='plasma', show=False)
plt.savefig("{0}/{1}_UMAP_IntestinalEnteroendocrineCell_marker_genes.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')

# Finding marker genes
# Rank genes using logistic regression. For instance, this has been suggested by Natranos et al. (2018). 
# The essential difference is that here, we use a multi-variate appraoch whereas conventional differential tests are uni-variate. Clark et al. (2014) 
sc.tl.rank_genes_groups(adata, 'louvain', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
plt.savefig("{0}/{1}_marker_gene_ranking.png".format(output_dir, bname) , bbox_inches='tight'); plt.close('all')


########################################################################
# Get the dataframe
adataDF = adata.to_df()

# Expression only subplot grid
nNeighbors  = [5, 12]
minDistance = [0.25,0.5]

import umap
selectedGenePathwayGenes = ['CHGA','CHGB','TAC1','TPH1','NEUROG3','SYP']
selectedGenePathwayGenes = markersDF['MarkerGenes'].tolist()
selectedGenePathwayGenes = ['HUMANMYC']
markerShapes = ["o", "X", "*", "v", "x", "d"]
colors = ['green', 'orangered']
for selectedGene in selectedGenePathwayGenes:
	try:
		# selectedGene = "Setd1b"
		# Convert anything greater than 0 to 1 .. this is used only for colors
		selectedGeneExp = np.array(adataDF[selectedGene].astype(int))
		selectedGeneExp = np.where(selectedGeneExp > 0, 'Yes', 'No')
		
		# Reset random state
		np.random.seed(2105)

		fig, ax = plt.subplots(len(nNeighbors), len(minDistance), sharex='col', sharey='row')
		for i in range(len(nNeighbors)):
			for j in range(len(minDistance)):
				n = nNeighbors[i]
				d = minDistance[j]
				np.random.seed(2105)
				reducer = umap.UMAP(n_neighbors = n, min_dist = d, random_state = 2105, transform_seed = 2105)
				embedding = reducer.fit_transform(np.array(adataDF))
				# Create the embedding dataframe
				embDF = pd.DataFrame(embedding, index=adataDF.index.values, columns=['UM1','UM2'])
				# Add cellType and expression information in the embedding DF
				embDF['Expressed'] = selectedGeneExp
				# plt.figure(figsize=(10,10))
				sns.scatterplot(x='UM1', y='UM2', palette = colors, hue='Expressed', hue_order =['Yes', 'No'], style  = 'Expressed', style_order =['Yes', 'No'], s = 5, data=embDF, edgecolor='none', alpha=0.75, ax=ax[i,j])
				ax[i,j].set_title("neighbors={0}, minDist={1}".format(n,d), fontsize=8)
				# Move the legend to an empty part of the plot
				ax[i,j].legend(loc='best', fontsize = 'x-small')
		plt.suptitle('UMAP projection of the {0} cells for {1} +/- '.format(projName, selectedGene), fontsize=12);
		# plt.show()
		# Output file
		output_file = "{0}/{1}_UMAP_marker_{2}.png".format(output_dir, bname, selectedGene)
		save_plot(output_file)
	except:
		print("\n- NOTE: {0} is not present in the dataset.".format(selectedGene))

