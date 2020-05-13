# pwd
cd /media/rad/SSD1/Genomes/STAR_References/MmHgMycIresCd2

# Add the humanMyc, IRES and human CD2 sequence to the mouse genome
cat /media/rad/SSD1/Genomes/STAR_References/mouse/GRCm38.primary_assembly.genome.fa hgMycIresCd2.fa > GRCm38_primary_assembly_hgMycIresCd2.fa
samtools faidx GRCm38_primary_assembly_hgMycIresCd2.fa

# Get the length of hgMycIresCd2.fa
samtools faidx hgMycIresCd2.fa
cut -f1-2 hgMycIresCd2.fa.fai
# hgMycIresCd2	2550

# Create combined human Myc IRES and CD2 GTF file (length = 2550)
# Names have to be the same
# Output shown in the docs section
cat hgMycIresCd2.gtf 

# Merge mouse and custom GTF 
cat /media/rad/SSD1/Genomes/STAR_References/mouse/gencode.vM22.annotation.gtf hgMycIresCd2.gtf > gencode_vM22_annotation_hgMycIresCd2.gtf

# Create STAR index
/home/rad/packages/STAR-2.7.0f/bin/Linux_x86_64/STAR --runThreadN 20 --runMode genomeGenerate --genomeDir vM22_hgMycIresCd2 --genomeFastaFiles GRCm38_primary_assembly_hgMycIresCd2.fa --sjdbGTFfile gencode_vM22_annotation_hgMycIresCd2.gtf


# at WS4
java -jar /home/rad/packages/picard_2.18.0/picard.jar \
  CreateSequenceDictionary \
  REFERENCE=GRCm38_primary_assembly_hgMycIresCd2.fa \
  OUTPUT=GRCm38_primary_assembly_hgMycIresCd2.dict

# refFlat
/home/rad/packages/ExportDropSeq/Drop-seq_tools-1.12/./ConvertToRefFlat \
  ANNOTATIONS_FILE=gencode_vM22_annotation_hgMycIresCd2.gtf \
  SEQUENCE_DICTIONARY=GRCm38_primary_assembly_hgMycIresCd2.dict \
  OUTPUT=GRCm38_primary_assembly_hgMycIresCd2.fa.refFlat

# Run cell ranger

# Copy H5 files to input folder
for t in bulk997 bulk1001 bulk1018 stomach1001;
do 
 rsync -arvP /media/rad/HDD2/temp_manec/${t}_mouse_MmHgMycIresCd2/outs/filtered_feature_bc_matrix.h5 /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/input/manec/tissue/${t}_mouse_filtered_feature_bc_matrix.h5
done

# Get the location of humanMyc, gap, ires and humanCD2 sequences in the hgMycIresCd2.fa file
# hgMycIresCd2.find(humanMyc) # The output is the location of first occurance of humanMyc in hgMycIresCd2
# In [8]: hgMycIresCd2.find(hmyc)
# Out[8]: 0

# In [9]: hgMycIresCd2.find(gap)
# Out[9]: 313

# In [10]: hgMycIresCd2.find(ires)
# Out[10]: 424

# In [11]: hgMycIresCd2.find(hcd2)
# Out[11]: 985

# Get the cell ID with hgMycIresCd2
for t in bulk997 bulk1001 bulk1018 stomach1001 bulk1079;
do 
  # Define sampledir
  sampledir="/media/rad/HDD2/temp_manec/${t}_mouse_MmHgMycIresCd2/outs"
  echo ${sampledir}
  # # Get filtered cell barcodes
  # cp ${sampledir}/filtered_feature_bc_matrix/barcodes.tsv.gz ${sampledir}
  # gunzip ${sampledir}/barcodes.tsv.gz

  # Get cell barcodes by:
  # - Subsample bam for hgMycIresCd2 region, 
  # - Filter for nM tag which 
  samtools view -b ${sampledir}/possorted_genome_bam.bam 'hgMycIresCd2:1-2550'| bamtools filter -tag nM:i:0 -in - | samtools view -h | LC_ALL=C grep -F -f ${sampledir}/barcodes.tsv | datamash transpose --no-strict -W | grep 'CB' | datamash transpose --no-strict -W | grep -Po 'CB:Z:[ACTG\-0-9]*'| sed 's/CB:Z://'| sort -u > ${sampledir}/${t}_hgMycIresCd2_cellIDs.txt

  samtools view -b ${sampledir}/possorted_genome_bam.bam 'hgMycIresCd2:1-313'| bamtools filter -tag nM:i:0 -in - | samtools view -h | LC_ALL=C grep -F -f ${sampledir}/barcodes.tsv | datamash transpose --no-strict -W | grep 'CB' | datamash transpose --no-strict -W | grep -Po 'CB:Z:[ACTG\-0-9]*'| sed 's/CB:Z://'| sort -u > ${sampledir}/${t}_hgMycIresCd2_humanMyc_cellIDs.txt

  # Filter for nM tag with mismatches more than 1
  samtools view ${sampledir}/possorted_genome_bam.bam '15:61985922-61989908'| egrep -vb "nM:i:0|nM:i:1" | LC_ALL=C grep -F -f ${sampledir}/barcodes.tsv | datamash transpose --no-strict -W | grep 'CB' | datamash transpose --no-strict -W | grep -Po 'CB:Z:[ACTG\-0-9]*'| sed 's/CB:Z://'| sort -u > ${sampledir}/${t}_hgMycIresCd2_humanMycMappedToMouseMyc_cellIDs.txt

  samtools view -b ${sampledir}/possorted_genome_bam.bam 'hgMycIresCd2:313-424'| bamtools filter -tag nM:i:0 -in - | samtools view -h | LC_ALL=C grep -F -f ${sampledir}/barcodes.tsv | datamash transpose --no-strict -W | grep 'CB' | datamash transpose --no-strict -W | grep -Po 'CB:Z:[ACTG\-0-9]*'| sed 's/CB:Z://'| sort -u > ${sampledir}/${t}_hgMycIresCd2_gap_cellIDs.txt

  samtools view -b ${sampledir}/possorted_genome_bam.bam 'hgMycIresCd2:424-985'| bamtools filter -tag nM:i:0 -in - | samtools view -h | LC_ALL=C grep -F -f ${sampledir}/barcodes.tsv | datamash transpose --no-strict -W | grep 'CB' | datamash transpose --no-strict -W | grep -Po 'CB:Z:[ACTG\-0-9]*'| sed 's/CB:Z://'| sort -u > ${sampledir}/${t}_hgMycIresCd2_ires_cellIDs.txt

  samtools view -b ${sampledir}/possorted_genome_bam.bam 'hgMycIresCd2:985-2550'| bamtools filter -tag nM:i:0 -in - | samtools view -h | LC_ALL=C grep -F -f ${sampledir}/barcodes.tsv | datamash transpose --no-strict -W | grep 'CB' | datamash transpose --no-strict -W | grep -Po 'CB:Z:[ACTG\-0-9]*'| sed 's/CB:Z://'| sort -u > ${sampledir}/${t}_hgMycIresCd2_humanCd2_cellIDs.txt
done

# 73 /media/rad/HDD2/temp_manec/bulk997_mouse_MmHgMycIresCd2/outs/bulk997_hgMycIresCd2_humanMycMappedToMouseMyc_cellIDs.txt
# 192 /media/rad/HDD2/temp_manec/bulk1001_mouse_MmHgMycIresCd2/outs/bulk1001_hgMycIresCd2_humanMycMappedToMouseMyc_cellIDs.txt
# 112 /media/rad/HDD2/temp_manec/bulk1018_mouse_MmHgMycIresCd2/outs/bulk1018_hgMycIresCd2_humanMycMappedToMouseMyc_cellIDs.txt
# 25 /media/rad/HDD2/temp_manec/stomach1001_mouse_MmHgMycIresCd2/outs/stomach1001_hgMycIresCd2_humanMycMappedToMouseMyc_cellIDs.txt

# Merge cell IDs into one file
cat /media/rad/HDD2/temp_manec/{bulk997,bulk1001,bulk1018,stomach1001,bulk1079}_mouse_MmHgMycIresCd2/outs/*_hgMycIresCd2_cellIDs.txt /media/rad/HDD2/temp_manec/{bulk997,bulk1001,bulk1018,stomach1001,bulk1079}_mouse_MmHgMycIresCd2/outs/*_hgMycIresCd2_humanMycMappedToMouseMyc_cellIDs.txt | sort -u > /media/rad/HDD2/temp_manec/hgMycIresCd2_cellIDs.txt

# 1832 /media/rad/HDD2/temp_manec/hgMycIresCd2_cellIDs.txt

for t in  humanMyc humanMycMappedToMouseMyc gap ires humanCd2;
do 
  cat /media/rad/HDD2/temp_manec/{bulk997,bulk1001,bulk1018,stomach1001,bulk1079}_mouse_MmHgMycIresCd2/outs/*_hgMycIresCd2_${t}_cellIDs.txt > /media/rad/HDD2/temp_manec/hgMycIresCd2_${t}_cellIDs.txt
  wc -l /media/rad/HDD2/temp_manec/hgMycIresCd2_${t}_cellIDs.txt
  for g in bulk997 bulk1001 bulk1018 stomach1001 bulk1079;
  do
    cat /media/rad/HDD2/temp_manec/${g}_mouse_MmHgMycIresCd2/outs/*_hgMycIresCd2_${t}_cellIDs.txt > /media/rad/HDD2/temp_manec/${g}_hgMycIresCd2_${t}_cellIDs.txt
    wc -l /media/rad/HDD2/temp_manec/${g}_hgMycIresCd2_${t}_cellIDs.txt
    
  done
done

# 431 /media/rad/HDD2/temp_manec/hgMycIresCd2_humanMyc_cellIDs.txt
# 709 /media/rad/HDD2/temp_manec/hgMycIresCd2_humanMycMappedToMouseMyc_cellIDs.txt
# 1597 /media/rad/HDD2/temp_manec/hgMycIresCd2_gap_cellIDs.txt
# 4108 /media/rad/HDD2/temp_manec/hgMycIresCd2_ires_cellIDs.txt
# 5407 /media/rad/HDD2/temp_manec/hgMycIresCd2_humanCd2_cellIDs.txt

################## DOCS #####################
# # cut -f1-2 GRCm38_primary_assembly_hgMycIresCd2.fa.fai
# 1	195471971
# 2	182113224
# 3	160039680
# 4	156508116
# 5	151834684
# 6	149736546
# 7	145441459
# 8	129401213
# 9	124595110
# 10	130694993
# 11	122082543
# 12	120129022
# 13	120421639
# 14	124902244
# 15	104043685
# 16	98207768
# 17	94987271
# 18	90702639
# 19	61431566
# X	171031299
# Y	91744698
# M	16299
# GL456210.1	169725
# GL456211.1	241735
# GL456212.1	153618
# GL456213.1	39340
# GL456216.1	66673
# GL456219.1	175968
# GL456221.1	206961
# GL456233.1	336933
# GL456239.1	40056
# GL456350.1	227966
# GL456354.1	195993
# GL456359.1	22974
# GL456360.1	31704
# GL456366.1	47073
# GL456367.1	42057
# GL456368.1	20208
# GL456370.1	26764
# GL456372.1	28664
# GL456378.1	31602
# GL456379.1	72385
# GL456381.1	25871
# GL456382.1	23158
# GL456383.1	38659
# GL456385.1	35240
# GL456387.1	24685
# GL456389.1	28772
# GL456390.1	24668
# GL456392.1	23629
# GL456393.1	55711
# GL456394.1	24323
# GL456396.1	21240
# JH584292.1	14945
# JH584293.1	207968
# JH584294.1	191905
# JH584295.1	1976
# JH584296.1	199368
# JH584297.1	205776
# JH584298.1	184189
# JH584299.1	953012
# JH584300.1	182347
# JH584301.1	259875
# JH584302.1	155838
# JH584303.1	158099
# JH584304.1	114452
# hgMycIresCd2	2550

# 2) cat hgMycIresCd2.gtf
# ┌──────────────┬────────┬────────────┬───┬──────┬───┬───┬───┬──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
# │ hgMycIresCd2 │ Custom │ gene       │ 1 │ 2550 │ . │ - │ . │ gene_id "hgMycIresCd2"; gene_type "hgMycIresCd2_inMouse_gene";  gene_name "hgMycIresCd2";                                                                                                                                                        │
# │ hgMycIresCd2 │ Custom │ transcript │ 1 │ 2550 │ . │ - │ . │ gene_id "hgMycIresCd2"; transcript_id "hgMycIresCd2"; gene_type "hgMycIresCd2_inMouse_gene"; gene_name "hgMycIresCd2"; transcript_type "hgMycIresCd2_inMouse_transcript"; transcript_name "hgMycIresCd2"; exon_number 1; exon_id "hgMycIresCd2"; │
# │ hgMycIresCd2 │ Custom │ exon       │ 1 │ 2550 │ . │ - │ . │ gene_id "hgMycIresCd2"; transcript_id "hgMycIresCd2"; gene_type "hgMycIresCd2_inMouse_gene"; gene_name "hgMycIresCd2"; transcript_type "hgMycIresCd2_inMouse_transcript"; transcript_name "hgMycIresCd2"; exon_number 1; exon_id "hgMycIresCd2"; │
# └──────────────┴────────┴────────────┴───┴──────┴───┴───┴───┴──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

