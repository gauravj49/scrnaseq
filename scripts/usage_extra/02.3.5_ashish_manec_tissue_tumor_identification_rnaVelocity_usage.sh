# pwd


# Get the annotation of spliced and unspliced genes
# Source: http://velocyto.org/velocyto.py/tutorial/cli.html#run10x-run-on-10x-chromium-samples
# Download repeatmasker file from https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Mouse&db=mm10&hgta_group=allTracks&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56694976-56714605&hgta_outputType=primaryTable&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf
rmaskGTF="/media/rad/HDD2/temp_manec/Ref_mouse_sanger_MmHgMycIresCd2/rmask/mm10_rmsk.gtf"
genesGTF="/media/rad/HDD2/temp_manec/Ref_mouse_sanger_MmHgMycIresCd2/genes/genes.gtf"

for t in bulk997 bulk1001 bulk1018 stomach1001;
do 
  echo "Processing: ${t}"
  sampleDir="/media/rad/HDD2/temp_manec/${t}_mouse_MmHgMycIresCd2"
  velocyto run10x -m ${rmaskGTF} ${sampleDir} ${genesGTF}
  echo ""
done

# Merging multiple samples/lanes in a single file
ipython

import loompy
files=["/media/rad/HDD2/temp_manec/bulk997_mouse_MmHgMycIresCd2/velocyto/bulk997_mouse_MmHgMycIresCd2.loom","/media/rad/HDD2/temp_manec/bulk1001_mouse_MmHgMycIresCd2/velocyto/bulk1001_mouse_MmHgMycIresCd2.loom","/media/rad/HDD2/temp_manec/bulk1018_mouse_MmHgMycIresCd2/velocyto/bulk1018_mouse_MmHgMycIresCd2.loom","/media/rad/HDD2/temp_manec/stomach1001_mouse_MmHgMycIresCd2/velocyto/stomach1001_mouse_MmHgMycIresCd2.loom"]
output_filename='/media/rad/HDD2/temp_manec/manec_tissues_merged_except1079_velocyto.loom'
loompy.combine(files, output_filename, key="Accession")

Cltr+D+D

# Run scVelo
ipython 

import scvelo as scv

# For beautified visualization you can change the matplotlib settings to our defaults with:
scv.settings.set_figure_params('scvelo')

# Read your data file (loom, h5ad, csv, …) using:
filename = '/media/rad/HDD2/temp_manec/manec_tissues_merged_except1079_velocyto.loom'
adata = scv.read(filename, cache=True)

# show proportions of spliced/unspliced abundances
scv.utils.show_proportions(adata)
adata

# Preprocess the data
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

# Remove duplicate cells
scv.pp.remove_duplicate_cells(adata)

# Calculate moments
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Compute velocity and velocity graph
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

# Plot results
scv.tl.umap(adata, random_state = 2105, n_components=3)
scv.pl.velocity_embedding_stream(adata, basis='umap')
scv.pl.velocity_embedding(adata, basis='umap', arrow_length=2, arrow_size=1.5, dpi=150)