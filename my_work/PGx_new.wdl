

version 1.0
#You need to declaration version information(version 1.0)
workflow PGx{
    input {
        File vcf_allsite_file
        File vcf_allsite_file_index
        #File reference_panel_dir 
        #File genetic_map_dir
        File bam_file
        File hla_file
        File bam_index
        
    }
    String sample_id = basename(vcf_allsite_file, ".allsite.vcf.gz")

    Array[String] chromosomes = ["chr1", "chr2","chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                                 "chr19", "chr20","chr21", "chr22","chrX"]
    # Array[Int] chromosomes = [20]

    # call allsite2toy{
    #     input:
    #         vcf_allsite_file = vcf_allsite_file,
    #         vcf_allsite_file_index = vcf_allsite_file_index
    # }
    
    scatter (chr in chromosomes) {
        
        call split_vcf{
            input:
                vcf_allsite_file  = vcf_allsite_file,
                vcf_allsite_file_index = vcf_allsite_file_index,
                chr = chr
        }

        call shapeit4_phase_task {
            input:
                vcf_file = split_vcf.vcf_file,
                vcf_file_index = split_vcf.vcf_file_index,
                reference_panel = split_vcf.reference_panel,
                reference_panel_index = split_vcf.reference_panel_index,
                chr_genetic_map = split_vcf.genetic_map,
                threads = 2,
                chr=chr
        }

    }

    call merge_phased_vcfs_task {
        input:
            
            phased_vcfs = shapeit4_phase_task.phased_vcf,
            vcf_file = split_vcf.vcf_file,
            unphased_vcfs_index = split_vcf.vcf_file_index
    }
    
    call cyriusTask{
        input:
            bam_file = bam_file,
            hla_file = hla_file,
            bam_index = bam_index
    }

    call preVcfPharmcatTask{
        input:
            vcf_file = merge_phased_vcfs_task.merged_phased_vcf,
            sample_id = sample_id
         
    }

    call pharmCatTask{
        input:
            file = preVcfPharmcatTask.pharmcat_vcf,
            outfile=cyriusTask.outcall_file
    }

    call PostPharmCat {

        # File match_report = "${filename}.match.json"
        # File phenotype_result = "${filename}.phenotype.json"
        # File report_file = "${filename}.report.json"
        # File report_tsv_result = "${filename}.report.tsv"
        
        input:
            report_json_file = pharmCatTask.report_file,
            report_tsv_file = pharmCatTask.report_tsv_result,
            
            outcall_file = cyriusTask.outcall_file,
            pre_pharmcat_vcf = preVcfPharmcatTask.pharmcat_vcf,
            merged_unphased_vcf = merge_phased_vcfs_task.merged_unphased_vcf,
            match_json_file = pharmCatTask.match_report,
            sample_id = sample_id
    }
  
  #String filename=basename(file,".vcf")
    output{
    #You need to define the workflow output to the “output” code block
        
        File match_report = pharmCatTask.match_report
        File phenotype_result = pharmCatTask.phenotype_result
        File report_tsv_result = pharmCatTask.report_tsv_result
        File report_file = pharmCatTask.report_file
        #File merged_phased_vcf = merge_phased_vcfs_task.merged_phased_vcf
        File pharmcat_vcf = preVcfPharmcatTask.pharmcat_vcf
        #File merged_unphased_vcf = merge_phased_vcfs_task.merged_unphased_vcf
        File outcall_file = cyriusTask.outcall_file
        
        File post_pharmcat_file = PostPharmCat.output_file
    }
  
 }

task split_vcf{
    input{
        File vcf_allsite_file
        File vcf_allsite_file_index
        String chr
        
    }

    command <<<
        ls /home/stereonote/map/
        source activate base
        conda activate bcftools
        bcftools view -o "~{chr}.allsite.vcf.gz" -O z -r "~{chr}" "~{vcf_allsite_file}"
        bcftools view -i 'INFO/AC > 0' ~{chr}.allsite.vcf.gz -Oz -o ~{chr}.vcf.gz
        bcftools index ~{chr}.vcf.gz
        cp /home/stereonote/map/~{chr}.b38.gmap.gz ~{chr}.b38.gmap.gz
        ref_dir="/zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/"
        ref_chr=""
        
    >>>
    output {
        File vcf_file = "~{chr}.vcf.gz"
        File vcf_file_index = "~{chr}.vcf.gz.csi"
        File genetic_map = "~{chr}.b38.gmap.gz"
        File reference_panel = "/zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_~{chr}.filtered.shapeit2-duohmm-phased.vcf.gz"    
        File reference_panel_index = "/zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_~{chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi"    
    }
    
    runtime {
    #The fixed parameter "docker_url" is used in wdl to specify the image address,Please copy the real url from "Image" and paste it here ,such as"stereonote/stereonote:latest". If you switch the area,please change the image url
    docker_url: "stereonote_hpc/longrui_c6c88994f0244878975c9042eeecb614_private:latest"
    #The fixed parameter " req_cpu" is used in wdl to apply for a task running cpu
    req_cpu: 1
    #The fixed parameter " req_memory" is used in wdl to specify the running memory
    req_memory: "1Gi"
    #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
  }
    
}

task shapeit4_phase_task {
    input {
        
        File vcf_file
        File vcf_file_index
        File reference_panel
        File reference_panel_index
        File chr_genetic_map
        Int threads = 1
        String chr
        
    }
    String vcf_file_name = basename(vcf_file, ".vcf.gz")

    command <<<
        
        echo ~{vcf_file}
        echo ~{chr_genetic_map}
        echo ~{threads},~{chr}
        echo ~{vcf_file_name}
        echo ~{reference_panel}
        /shapeit4/shapeit4/bin/shapeit4.2 \
            --input ~{vcf_file} \
            --reference ~{reference_panel} \
            --map ~{chr_genetic_map} \
            --output ~{vcf_file_name}.phased.vcf \
            --thread ~{threads} \
            --sequencing \
            --region ~{chr} 
    >>>

    output {
        File phased_vcf = "~{vcf_file_name}.phased.vcf"
        
    }

    runtime {
        docker_url: "stereonote_hpc_external/longrui_9a9c6fdbbbe34ecbb3d83a24c507b316_private:latest"
        req_cpu: 1
        #The fixed parameter " req_memory" is used in wdl to specify the running memory
        req_memory: "2.5Gi"
    }
}
task merge_phased_vcfs_task {
    input {
        #File phase_dir
        Array[File] phased_vcfs
        Array[File] vcf_file
        Array[File] unphased_vcfs_index
    }
    
    #String phased_vcf_files = sep=" " phased_vcfs
    command {
        
        # bcftools concat -o merged_phased.vcf -O v $phased_vcf_files
        bcftools concat -o merged_phased.vcf -O v ${sep=' ' phased_vcfs}
        bcftools concat -o merged_unphased.vcf -O v  ${sep=' ' vcf_file}
    }
    output {
        File merged_phased_vcf = "merged_phased.vcf"
        File merged_unphased_vcf = "merged_unphased.vcf"
        
    }
    runtime {
        docker_url: "public-library/shuliping_b2db8c5c625743c79ca047f52ec29070_public:latest"
        req_cpu: 1
        req_memory: "1Gi"
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
    java -jar /home/stereonote/pharmcat-3.0.1-all.jar -vcf ~{file} -po ~{outfile}  -reporterJson --output-dir .  -reporterCallsOnlyTsv
    
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
    #pharmcat.example.report.html
  }
}
task cyriusTask{
    input {
        File bam_file
        File hla_file
        File bam_index
    }
    command {
        export PATH=/home/stereonote/script:$PATH
        echo $PATH
        echo "hello,world "
        echo ~{bam_file}>bam_path.txt
        
        # docker run --rm \
        cyrius --manifest bam_path.txt \
                       --genome 38 \
                       --prefix cyp2d6 \
                       --outDir . \
                       --threads 4
        cyp2d6_results_handler.py --input_file cyp2d6.tsv
        hla_results_handler.py --input_file ~{hla_file}
        cat cyp2d6_result.txt>outcall.tsv
        echo >> outcall.tsv
        cat hla_result.txt>>outcall.tsv
    }
    runtime {
        #The fixed parameter "docker_url" is used in wdl to specify the image address,Please copy the real url from "Image" and paste it here ,such as"stereonote/stereonote:latest". If you switch the area,please change the image url
        docker_url: "stereonote_hpc/longrui_bf9939ec3d4341378bfb30423fc5b1ef_private:latest"
        
        
        #The fixed parameter " req_cpu" is used in wdl to apply for a task running cpu
        req_cpu: 4
        #The fixed parameter " req_memory" is used in wdl to specify the running memory
        req_memory: "2Gi"
        #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
        
    }
    output {
        File outcall_file = "outcall.tsv"

    }
}

task preVcfPharmcatTask{
    input {
        File vcf_file
        String sample_id
      
    }
    
    command <<<
        
        echo ~{vcf_file}
        # copy and paste the command to the docker image
        mkdir -p results/pharmcat_ready/
        source activate /home/stereonote/miniconda3
        /home/stereonote/miniconda3/bin/python /home/stereonote/preprocessor/pharmcat_vcf_preprocessor.py \
        -vcf ~{vcf_file} \
        -refFna /home/stereonote/ref_files/reference.fna.bgz \
        -refVcf /home/stereonote/ref_files/pharmcat_positions_3.0.1.vcf.bgz \
        -o results/pharmcat_ready/ \
        -bf ~{sample_id}
        bgzip -d -c results/pharmcat_ready/~{sample_id}.preprocessed.vcf.bgz > ~{sample_id}.preprocessed.vcf
    >>>

    output {
        File pharmcat_vcf_bgz = glob("results/pharmcat_ready/*.vcf.bgz")[0]
        File pharmcat_vcf = "~{sample_id}.preprocessed.vcf"
    }

    runtime {
        docker_url: "stereonote_hpc/longrui_f578730aae604484ac1296f23e069cd3_private:latest"
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
        File match_json_file
        String sample_id
    }

    command <<<
        source activate base
        conda activate bcftools
        java -Xmx2g -jar /opt/conda/envs/bcftools/share/snpeff-5.2-1/snpEff.jar -v GRCh38.105  ~{pre_pharmcat_vcf} >annotated.vcf
        bgzip annotated.vcf
        bcftools index annotated.vcf.gz
        bgzip ~{merged_unphased_vcf}
        bcftools index ~{merged_unphased_vcf}.gz
        mkdir -p isec
        bcftools isec  -p isec -Oz -W ~{merged_unphased_vcf}.gz  annotated.vcf.gz
        bcftools view  isec/0002.vcf.gz| grep -v '^##' > intersect.csv
        bcftools view  annotated.vcf.gz| grep -v '^##' > annotated.csv
        sed -i 's/#CHROM/CHROM/g' intersect.csv
        sed -i 's/#CHROM/CHROM/g' annotated.csv

        # bgzip ~{pre_pharmcat_vcf}
        # bcftools index ~{pre_pharmcat_vcf}.gz
        # bgzip ~{merged_unphased_vcf}
        # bcftools index ~{merged_unphased_vcf}.gz
        # bcftools isec  -p isec -Oz -W ~{merged_unphased_vcf}.gz  ~{pre_pharmcat_vcf}.gz
        # bcftools view  isec/0002.vcf.gz| grep -v '^##' > intersect.csv
        # bcftools view  ~{pre_pharmcat_vcf}.gz| grep -v '^##' > pre_pharmcat.csv
        # sed -i 's/#CHROM/CHROM/g' intersect.csv
        # sed -i 's/#CHROM/CHROM/g' pre_pharmcat.csv
        python /home/stereonote/script/post_pharmcat.py \
            -report_json_file ~{report_json_file} \
            -result_tsv_file ~{report_tsv_file} \
            -outcall_file ~{outcall_file} \
            -intersect_file intersect.csv \
            -prepharmcat_file annotated.csv \
            -sample_id ~{sample_id} \
    >>>

    output {
        File output_file = "~{sample_id}.pgx.tsv"
    }

    runtime {
        docker_url: "stereonote_hpc/longrui_ccc92f50b0b14dad83dc198bcb4646ea_private:latest"
        req_cpu: 1
        req_memory: "3Gi"
  }
}