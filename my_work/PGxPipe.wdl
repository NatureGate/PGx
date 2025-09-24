

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
                                 "chr19", "chr20","chr21", "chr22"]
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
                threads = 8,
                chr=chr
        }

    }

    

    call merge_phased_vcfs_task {
        input:
            
            phased_vcfs = shapeit4_phase_task.phased_vcf,
            unphased_vcfs = split_vcf.vcf_file,
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
            vcf_file = merge_phased_vcfs_task.merged_phased_vcf
            
            
    }

    call pharmCatTask{
        input:
            file = preVcfPharmcatTask.pharmcat_vcf,
            outfile=cyriusTask.outcall_file
    }
  
  #String filename=basename(file,".vcf")
    output{
    #You need to define the workflow output to the “output” code block
        
        File match_report = pharmCatTask.match_report
        File phenotype_result = pharmCatTask.phenotype_result
        File report_file = pharmCatTask.report_file
        File merged_phased_vcf = merge_phased_vcfs_task.merged_phased_vcf
        File pharmcat_vcf = preVcfPharmcatTask.pharmcat_vcf
        File merged_unphased_vcf = merge_phased_vcfs_task.merged_unphased_vcf
        
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
    req_cpu: 4
    #The fixed parameter " req_memory" is used in wdl to specify the running memory
    req_memory: "8Gi"
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
        Int threads = 8
        String chr
        
    }
    String vcf_file_name = basename(vcf_file, ".vcf.gz")

    command <<<
        
        echo ~{vcf_file}
        echo ~{chr_genetic_map}
        echo ~{threads},~{chr}
        echo ~{vcf_file_name}
        echo ~{reference_panel}
        shapeit4 \
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
        docker_url: "stereonote_hpc_external/longrui_3895463e051d4de59b736f5b66870bea_private:latest"
        req_cpu: 8
        #The fixed parameter " req_memory" is used in wdl to specify the running memory
        req_memory: "4Gi"
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
        bcftools concat -o merged_unphased.vcf -O v ${sep=' ' vcf_file}
    }
    output {
        File merged_phased_vcf = "merged_phased.vcf"
        File merged_unphased_vcf = "merged_unphased.vcf"
    }
    runtime {
        docker_url: "public-library/shuliping_b2db8c5c625743c79ca047f52ec29070_public:latest"
        req_cpu: 1
        req_memory: "4Gi"
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
    java -jar /home/stereonote/pharmcat-3.0.0-all.jar -vcf ~{file} -po ~{outfile}  -reporterJson --output-dir .  
    ls -l
    pwd
  }
  runtime {
    #The fixed parameter "docker_url" is used in wdl to specify the image address,Please copy the real url from "Image" and paste it here ,such as"stereonote/stereonote:latest". If you switch the area,please change the image url
    docker_url: "stereonote_hpc/longrui_d8440a5f4e2245c29c01e170f20c2682_private:latest"
    #The fixed parameter " req_cpu" is used in wdl to apply for a task running cpu
    req_cpu: 2
    #The fixed parameter " req_memory" is used in wdl to specify the running memory
    req_memory: "4Gi"
    #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
  }
  output {
    File match_report = "${filename}.match.json"
    File phenotype_result = "${filename}.phenotype.json"
    File report_file = "${filename}.report.json"
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
                       --threads 8
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
        req_cpu: 8
        #The fixed parameter " req_memory" is used in wdl to specify the running memory
        req_memory: "10Gi"
        #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
        
    }
    output {
        File outcall_file = "outcall.tsv"

    }
}

task allsite2toy{
    input{
        File vcf_allsite_file
        File vcf_allsite_file_index
    }

    command <<<
        source activate base
        conda activate bcftools
        echo ~{vcf_allsite_file}
        echo ~{vcf_allsite_file_index}
        
        sh /home/stereonote/script/allsite2small.sh ~{vcf_allsite_file}
        
    >>>
    output {
        File toy_vcf_file = "final_merged/merged_all.vcf.gz.sorted.gz"
        File toy_vcf_file_index = "final_merged/merged_all.vcf.gz.sorted.gz.csi"
    }
    
    runtime {
    #The fixed parameter "docker_url" is used in wdl to specify the image address,Please copy the real url from "Image" and paste it here ,such as"stereonote/stereonote:latest". If you switch the area,please change the image url
    docker_url: "stereonote_hpc/longrui_1f89c93319094898ae90623727f1d0d3_private:latest"
    #The fixed parameter " req_cpu" is used in wdl to apply for a task running cpu
    req_cpu: 4
    #The fixed parameter " req_memory" is used in wdl to specify the running memory
    req_memory: "8Gi"
    #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
  }
    
}

task preVcfPharmcatTask{
    input {
        File vcf_file
      
    }
    
    command <<<
        
        echo ~{vcf_file}
        # copy and paste the command to the docker image
        mkdir -p results/pharmcat_ready/
        python3 /pharmcat/pharmcat_vcf_preprocessor.py \
        -vcf ~{vcf_file} \
        -refFna /pharmcat/reference.fna.bgz \
        -refVcf /pharmcat/pharmcat_positions.vcf.bgz \
        -o results/pharmcat_ready/
    >>>

    output {
        File pharmcat_vcf = glob("results/pharmcat_ready/*.vcf")[0]
    }

    runtime {
        docker_url: "stereonote_hpc_external/longrui_8725a03a53424685888b56480b3f7f08_private:latest"
        req_cpu: 1
        #The fixed parameter " req_memory" is used in wdl to specify the running memory
        req_memory: "2Gi"
    }
}


