# Run cell ranger with custom transcriptome

cellranger count --id=bulk1079_mouse_MmHgMycIresCd2 --transcriptome=/media/rad/HDD2/temp_manec/Ref_mouse_sanger_MmHgMycIresCd2 --fastqs=/media/rad/HDD2/temp_manec/10x_batch2/fastqs/HC223BGXC/10x_batch2_1079_bulk/ --sample=10x_batch2_1079_bulk --localcores=24 --expect-cells=10000

cellranger count --id=bulk1001_mouse_MmHgMycIresCd2 --transcriptome=/media/rad/HDD2/temp_manec/Ref_mouse_sanger_MmHgMycIresCd2 --fastqs=/media/rad/HDD2/temp_manec/10x_batch2/fastqs/HC223BGXC/10x_batch2_1001_bulk/ --sample=10x_batch2_1001_bulk --localcores=24 --expect-cells=7000

cellranger count --id=bulk1018_mouse_MmHgMycIresCd2 --transcriptome=/media/rad/HDD2/temp_manec/Ref_mouse_sanger_MmHgMycIresCd2 --fastqs=/media/rad/HDD2/temp_manec/10x_batch2/fastqs/HC223BGXC/10x_batch2_1018_bulk/ --sample=10x_batch2_1018_bulk --localcores=24 --expect-cells=7000

cellranger count --id=bulk997_mouse_MmHgMycIresCd2 --transcriptome=/media/rad/HDD2/temp_manec/Ref_mouse_sanger_MmHgMycIresCd2 --fastqs=/media/rad/HDD2/temp_manec/10x_batch2/fastqs/HC223BGXC/10x_batch2_997_bulk/ --sample=10x_batch2_997_bulk --localcores=24 --expect-cells=7000

cellranger count --id=stomach1001_mouse_MmHgMycIresCd2 --transcriptome=/media/rad/HDD2/temp_manec/Ref_mouse_sanger_MmHgMycIresCd2 --fastqs=/media/rad/HDD2/temp_manec/10x_batch2/fastqs/HC223BGXC/10x_batch2_1001_stomach/ --sample=10x_batch2_1001_stomach --localcores=24 --expect-cells=7000
