version 1.0
#You need to declaration version information(version 1.0)


workflow PGx{
    input {
        File outfile
        File merged_phased_vcf 
        File merged_unphased_vcf 
        File merged_phased_vcf_gz 
        File merged_phased_vcf_gz_index 
        File merged_unphased_vcf_gz 
        File merged_unphased_vcf_gz_index 
        String SampleID
        
        
        
    }
    

    call preVcfPharmcatTask{
        input:
            sample_id = SampleID,
            outfile = outfile,
            merged_phased_vcf = merged_phased_vcf,
            merged_unphased_vcf = merged_unphased_vcf,
            merged_phased_vcf_gz = merged_phased_vcf_gz,
            merged_phased_vcf_gz_index = merged_phased_vcf_gz_index,
            merged_unphased_vcf_gz = merged_unphased_vcf_gz,
            merged_unphased_vcf_gz_index = merged_unphased_vcf_gz_index

         
    }

    call pharmCatTask{
        input:
            file = preVcfPharmcatTask.pharmcat_vcf,
            outfile = preVcfPharmcatTask.update_outcall
    }

    call PostPharmCat {
        input:
            report_json_file = pharmCatTask.report_file,
            report_tsv_file = pharmCatTask.report_tsv_result,
            outcall_file = preVcfPharmcatTask.update_outcall,
            pre_pharmcat_vcf = preVcfPharmcatTask.pharmcat_vcf,
            merged_unphased_vcf = merged_unphased_vcf,
            merged_phased_vcf_gz = merged_phased_vcf_gz,
            merged_phased_vcf_gz_index = merged_phased_vcf_gz_index,
            match_json_file = pharmCatTask.match_report,
            sample_id = SampleID
    }
  
  #String filename=basename(file,".vcf")
    output{
        
        File match_report = pharmCatTask.match_report
        File phenotype_result = pharmCatTask.phenotype_result
        File report_tsv_result = pharmCatTask.report_tsv_result
        File report_file = pharmCatTask.report_file
        #File merged_phased_vcf = merge_phased_vcfs_task.merged_phased_vcf
        File pharmcat_vcf = preVcfPharmcatTask.pharmcat_vcf
        #File merged_unphased_vcf = merge_phased_vcfs_task.merged_unphased_vcf      
        File unknown_genotype_file = PostPharmCat.unknown_genotype_file
        File other_file = PostPharmCat.other_file
        File post_pharmcat_file = PostPharmCat.output_file
    }
 }

task pharmCatTask{
  input {
    
    File file
    String filename=basename(file,".preprocessed.vcf")
    File outfile
  }
  command {
  
    echo "hello,world "
    echo ${file},${filename}
    # docker run --rm \
    export JAVA_HOME=/home/stereonote/jdk17
    export CLASSPATH=/home/stereonote:/home/stereonote/jdk17/lib
    export PATH=/home/stereonote:/home/stereonote/jdk17/bin:$PATH
    source activate base
    # echo -e "HLA-A\t*31:01/*32:01\nHLA-B\t*15:02/*57:01\nCYP2D6\t*1/*3" > outcall.tsv
    java -jar /home/stereonote/pharmcat-3.0.1-all.jar -vcf ~{file}  -po ~{outfile} -reporterJson --output-dir .  -reporterCallsOnlyTsv

  }
  runtime {
    #The fixed parameter "docker_url" is used in wdl to specify the image address,Please copy the real url from "Image" and paste it here ,such as"stereonote/stereonote:latest". If you switch the area,please change the image url
    docker_url: "stereonote_hpc/longrui_f578730aae604484ac1296f23e069cd3_private:latest"
    #The fixed parameter " req_cpu" is used in wdl to apply for a task running cpu
    req_cpu: 1
    #The fixed parameter " req_memory" is used in wdl to specify the running memory
    req_memory: "2Gi"
    #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
  }
  output {
    File match_report = "${filename}.match.json"
    File phenotype_result = "${filename}.phenotype.json"
    File report_file = "${filename}.report.json"
    File report_tsv_result = "${filename}.report.tsv"
    # File outcall_file = "outcall.tsv"
    #pharmcat.example.report.html
  }
}


task preVcfPharmcatTask{
    input {
        
        String sample_id
        File outfile
        File merged_phased_vcf 
        File merged_unphased_vcf 
        File merged_phased_vcf_gz 
        File merged_phased_vcf_gz_index 
        File merged_unphased_vcf_gz 
        File merged_unphased_vcf_gz_index
        
      
    }
    
    String mt_rnr1_type = "rs267606618,rs267606619,rs267606617,rs28358569"
    
    command <<<
        source activate /home/stereonote/miniconda3
        bcftools filter -e 'FILTER="LowQual"' ~{merged_unphased_vcf_gz} -o filtered_merged_unphased.vcf.gz
        bcftools index filtered_merged_unphased.vcf.gz
        bcftools annotate -a ~{merged_phased_vcf_gz} -c FORMAT/GT -Ov -o annotated.vcf filtered_merged_unphased.vcf.gz
        bcftools view -i "ID == 'rs267606618' || ID == 'rs267606619' || ID == 'rs267606617' || ID == 'rs28358569'" ~{merged_unphased_vcf_gz}|grep -v '^#'|awk -F'\t' '{print $3}'|head -n1>MTRNR1.txt
        cp ~{outfile}   ~{sample_id}.outcall.tsv
        python /home/stereonote/script/check_MTRNR1.py MTRNR1.txt ~{sample_id}.outcall.tsv
        # copy and paste the command to the docker image
        mkdir -p results/pharmcat_ready/
        source activate /home/stereonote/miniconda3
        /home/stereonote/miniconda3/bin/python /home/stereonote/preprocessor/pharmcat_vcf_preprocessor.py \
        -vcf annotated.vcf \
        -refFna /home/stereonote/ref_files/reference.fna.bgz \
        -refVcf /home/stereonote/ref_files/pharmcat_positions_3.0.1.vcf.bgz \
        -o results/pharmcat_ready/ \
        -bf ~{sample_id}
        bgzip -d -c results/pharmcat_ready/~{sample_id}.preprocessed.vcf.bgz > ~{sample_id}.preprocessed.vcf
    >>>

    output {
        File pharmcat_vcf_bgz = glob("results/pharmcat_ready/*.vcf.bgz")[0]
        File pharmcat_vcf = "~{sample_id}.preprocessed.vcf"
        File update_outcall = "~{sample_id}.outcall.tsv"
    }

    runtime {
        docker_url: "stereonote_hpc/longrui_f425a0a82aed451bb21bc95f35021e07_private:latest"
        req_cpu: 1
        #The fixed parameter " req_memory" is used in wdl to specify the running memory
        req_memory: "1Gi"
    }
}

task PostPharmCat {
    input {
        File report_json_file
        File report_tsv_file
        File outcall_file
        File pre_pharmcat_vcf
        File merged_unphased_vcf
        File merged_phased_vcf_gz
        File merged_phased_vcf_gz_index
        File match_json_file
        String sample_id
    }

    command <<<
        source activate base
        conda activate bcftools
        java -Xmx2g -jar /opt/conda/envs/bcftools/share/snpeff-5.2-1/snpEff.jar -v GRCh38.105  ~{pre_pharmcat_vcf} >annotated.vcf
        bgzip annotated.vcf
        bgzip -c ~{merged_unphased_vcf}>merged_unphased.vcf.gz
        bcftools index annotated.vcf.gz
        bcftools index merged_unphased.vcf.gz
        
        mkdir -p isec
        bcftools isec  -p isec -Oz -W merged_unphased.vcf.gz  annotated.vcf.gz
        bcftools view  isec/0002.vcf.gz| grep -v '^##' > intersect.csv
        bcftools view  annotated.vcf.gz| grep -v '^##' > annotated.csv
        sed -i 's/#CHROM/CHROM/g' intersect.csv
        sed -i 's/#CHROM/CHROM/g' annotated.csv
        bcftools query -f '%ID\n' ~{merged_phased_vcf_gz} | grep -v '^.$' > valid_ids.txt
        
        python /home/stereonote/script/post_pharmcat_new.py \
            -report_json_file ~{report_json_file} \
            -result_tsv_file ~{report_tsv_file} \
            -outcall_file ~{outcall_file} \
            -intersect_file intersect.csv \
            -prepharmcat_file annotated.csv \
            -sample_id ~{sample_id} 
        # python /home/stereonote/script/var_drug.py --valid_ids_path valid_ids.txt --sample_id ~{sample_id}
    >>>

    output {
        File output_file = "~{sample_id}.pgx.tsv"
        File unknown_genotype_file = "~{sample_id}_unknown.pgx.tsv"
        File other_file = "~{sample_id}_other.pgx.tsv"
    }

    runtime {
        docker_url: "stereonote_hpc/longrui_674c466902cf4e009b945944ef4c4b16_private:latest"
        req_cpu: 1
        req_memory: "3Gi"
  }
}