#!/bin/bash
# This script is used to phase a VCF file using ShapeIt4.2
    # command = [
    #     'bash', script_path,
    #     chr,
    #     vcf_file,
        
    # ]
# 添加上面command中的参数

chr=$1
vcf_file=$2
vcf_file_index=${vcf_file}.csi
chr_genetic_map="/home/stereonote/map/${chr}.b38.gmap.gz"
reference_panel="/zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
reference_panel_index="/zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi"
echo "vcf_file: $vcf_file"
echo "vcf_file_index: $vcf_file_index"
echo "chr_genetic_map: $chr_genetic_map"
echo "reference_panel: $reference_panel"
echo "reference_panel_index: $reference_panel_index"
#map_dir=$2
threads=1

# reference_panel="/zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_~{chr}.filtered.shapeit2-duohmm-phased.vcf.gz"    
# reference_panel_index="/zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_~{chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi" 
/shapeit4/shapeit4/bin/shapeit4.2 \
            --input ${vcf_file} \
            --reference ${reference_panel} \
            --map ${chr_genetic_map} \
            --output ${chr}.phased.vcf \
            --thread ${threads} \
            --sequencing \
            --region ${chr} 