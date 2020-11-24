# PWD
cd /home/rad/users/gaurav/projects/seqAnalysis/scrnaseq

# Run cellranger with --force-cells option
cellranger count --id=10x_220819_1018_bulk --fastqs=/media/rad/SSD1/temp_manec/fastq/H7YYWBGXC/10x_220819_1018_bulk/ --sample=10x_220819_1018_bulk --output-dir=/media/rad/SSD1/temp_manec/output --expect-cells=10000 --transcriptome=/home/rad/users/ashish/refdata-cellranger-GRCh38-and-mm10-3.1.0 --localcores=4