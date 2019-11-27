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

