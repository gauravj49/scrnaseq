# pwd
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

projdir="/media/rad/HDD2/temp_manec/bulk1018_mouse_sangerhMYC_IRES/outs"
cd ${projdir}

samtools view -h ${projdir}/bulk1018_Myc.zeroEditDistance.bam | LC_ALL=C grep -F -f barcodes.tsv | datamash transpose --no-strict -W | grep 'CB' | datamash transpose --no-strict -W | sort -u > ${projdir}/bulk1018_Myc_human_cellIDs.txt

samtools view -h ${projdir}/bulk1018_IRESreads.bam | LC_ALL=C grep -F -f barcodes.tsv | datamash transpose --no-strict -W | grep 'CB' | datamash transpose --no-strict -W | sort -u > ${projdir}/bulk1018_IRESreads_cellIDs.txt



















# # Get the cell barcodes for every sample containing human myc
# subsetBamDir="/home/rad/packages/cellranger-3.0.2/subset-bam/subset-bam-1.0-x86_64-linux"
# ${subsetBamDir}/subset-bam --bam bulk1018_Myc.zeroEditDistance.bam --cell-barcodes filter.txt --out-bam test.bam


# samtools view -h ${projdir}/bulk1018_Myc.zeroEditDistance.bam | grep -w CB | cut -f26 | sed 's/CB\:Z\://g' | sort -u > ${projdir}/bulk1018_Myc_human_cellIDs.txt


