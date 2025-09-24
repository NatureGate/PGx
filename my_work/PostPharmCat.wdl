version 1.0
workflow{
    input {
        File report_json_file
        File report_tsv_file
        File outcall_file
        File pre_pharmcat_vcf
        File merged_unphased_vcf
        File match_json_file
        String sample_id

    }

    call PostPharmCat {
        input:
            report_json_file = report_json_file,
            report_tsv_file = report_tsv_file,
            outcall_file = outcall_file,
            pre_pharmcat_vcf = pre_pharmcat_vcf,
            merged_unphased_vcf = merged_unphased_vcf,
            match_json_file = match_json_file,
            sample_id = sample_id
            
    }

    output {
        File output_file = PostPharmCat.output_file
    }
}

task PostPharmCat {
    input {
        File report_json_file
        File report_tsv_file
        File outcall_file
        File pre_pharmcat_vcf
        File merged_unphased_vcf
        File match_json_file
        String sample_id
    }

    command <<<
        bgzip ~{pre_pharmcat_vcf}
        bcftools index ~{pre_pharmcat_vcf}.gz
        bgzip ~{merged_unphased_vcf}
        bcftools index ~{merged_unphased_vcf}.gz
        mkdir -p isec
        bcftools isec  -p isec -Oz -W ~{merged_unphased_vcf} ~{pre_pharmcat_vcf}.gz
        bcftools view  isec/0002.vcf.gz| grep -v '^##' > intersect.csv
        bcftools view  ~{pre_pharmcat_vcf}.gz| grep -v '^##' > pre_pharmcat.csv
        sed -i 's/#CHROM/CHROM/g' intersect.csv
        sed -i 's/#CHROM/CHROM/g' pre_pharmcat.csv
        source activate base
        conda activate bcftools
        python /home/stereonote/script/post_pharmcat.py \
            -report_json_file ~{report_json_file} \
            -result_tsv_file ~{report_tsv_file} \
            -outcall_file ~{outcall_file} \
            -intersect_file intersect.csv \
            -prepharmcat_file pre_pharmcat.csv
    >>>

    output {
        File output_file = "drug_gene_table.csv"
    }

    runtime {
        docker_url: "stereonote_hpc/longrui_ba500a9365b549dabcf990b8bc98affe_private:latest"
        req_cpu: 1
        req_memory: "1Gi"
  }
}