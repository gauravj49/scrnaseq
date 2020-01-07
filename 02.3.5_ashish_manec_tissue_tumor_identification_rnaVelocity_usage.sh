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
