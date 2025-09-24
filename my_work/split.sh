#!/bin/bash
source activate base
conda activate bcftools
chr=$1
vcf_file=$2
vcf_file_index=$3
echo "vcf_file: $vcf_file"
echo "vcf_file_index: $vcf_file_index"
#map_dir=$2
threads=1


bcftools view -o "${chr}.vcf.gz" -O z -r "${chr}" "${vcf_file}" 
bcftools index ${chr}.vcf.gz
cp /home/stereonote/map/${chr}.b38.gmap.gz ${chr}.b38.gmap.gz
# reference_panel="/zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
# reference_panel_index="/zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi"
