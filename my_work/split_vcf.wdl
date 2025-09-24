

version 1.0
#You need to declaration version information(version 1.0)
workflow PGx{
    input {
        File vcf_allsite_file
        File vcf_allsite_file_index
        File reference_panel_dir 
        File genetic_map_dir
        File bam_file
        String name=longrui

    }

    Array[String] chromosomes = ["chr1", "chr2","chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                                 "chr19", "chr20","chr21", "chr22", "chrX"]
    # Array[Int] chromosomes = [20]

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
                reference_panel_dir = reference_panel_dir,
                genetic_map_dir = genetic_map_dir,
                threads = 16,
                chr=chr
        }

    }
    call merge_phased_vcfs_task {
        input:
            
            phased_vcfs = shapeit4_phase_task.phased_vcf
    }
    # call shapeit4_phase {
    #     input:
    #         vcf_file = vcf_file_dir,
    #         vcf_file_index = vcf_file_index_dir,
    #         reference_panel = reference_panel_dir,
    #         reference_panel_index = reference_panel_index_dir,
    #         genetic_map = genetic_map_dir,
    #         threads = threads
    # }
     call pharmCatTask{
         input:
         file = merge_phased_vcfs_task.merged_phased_vcf
     }
  
  #String filename=basename(file,".vcf")
    output{
    #You need to define the workflow output to the “output” code block
        
        File match_report = pharmCatTask.match_report
        File phenotype_result = pharmCatTask.phenotype_result
        File report_file = pharmCatTask.report_file
        
    }
  
 }

task split_vcf{
    input{
        File vcf_allsite_file="/home/stereonote/test/merged_filtered_vcfs/sorted.vcf.gz"
        File vcf_allsite_file_index="/home/stereonote/test/merged_filtered_vcfs/sorted.vcf.gz.csi"
        String chr
        
    }

    command <<<
        source activate base
        conda activate bcftools
        echo ~{vcf_allsite_file}
        echo ~{vcf_allsite_file_index}
        echo ~{chr}
        
        bcftools view -o "~{chr}.vcf.gz" -O z -r "~{chr}" "~{vcf_allsite_file}"
        bcftools view -i 'INFO/AC > 0' ~{chr}.allsite.vcf.gz -Oz -o ~{chr}.vcf.gz
        bcftools index ~{chr}.vcf.gz
        
    >>>
    output {
        File vcf_file = "~{chr}.vcf.gz"
        File vcf_file_index = "~{chr}.vcf.gz.tbi"
    }
    
    runtime {
    #The fixed parameter "docker_url" is used in wdl to specify the image address,Please copy the real url from "Image" and paste it here ,such as"stereonote/stereonote:latest". If you switch the area,please change the image url
    docker_url: "stereonote_hpc/longrui_c6c88994f0244878975c9042eeecb614_private:latest"
    #The fixed parameter " req_cpu" is used in wdl to apply for a task running cpu
    req_cpu: 4
    #The fixed parameter " req_memory" is used in wdl to specify the running memory
    req_memory: "64Gi"
    #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
  }
    
}
task pharmCatTask{
  input {
    
    File file
    String filename=basename(file,".vcf")
  }
  command {
    echo "hello,world "
    echo ${file},${filename}
    # docker run --rm \
    pharmcat_pipeline --output-dir pharmcat ${file}
    ls -l
    pwd
  }
  runtime {
    #The fixed parameter "docker_url" is used in wdl to specify the image address,Please copy the real url from "Image" and paste it here ,such as"stereonote/stereonote:latest". If you switch the area,please change the image url
    docker_url: "public-library/longrui_4605c9b9e91b4812a4a9f696ec47ca7a_public:latest"
    #The fixed parameter " req_cpu" is used in wdl to apply for a task running cpu
    req_cpu: 2
    #The fixed parameter " req_memory" is used in wdl to specify the running memory
    req_memory: "4Gi"
    #The fixed parameter "docker url" is used in wdl to specify the image address,Please copy the image url from "Image" and paste it here
  }
  output {
    File match_report = "pharmcat/${filename}.match.json"
    File phenotype_result = "pharmcat/${filename}.phenotype.json"
    File report_file = "pharmcat/${filename}.report.html"
    #pharmcat.example.report.html
  }
}

task shapeit4_phase_task {
    input {
        
        File vcf_file
        File vcf_file_index
        File reference_panel_dir= "/zfsms3/pub/database/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/"
        File genetic_map_dir="/home/stereonote/map"
        Int threads = 8
        String chr
        
    }
    String vcf_file_name = basename(vcf_file, ".vcf.gz")

    command <<<
        ls -l
        reference_panel=$(ls ~{reference_panel_dir}|grep "_~{chr}.filtered.shapeit2-duohmm-phased.vcf.gz$")
        genetic_map=$(ls ~{genetic_map_dir}|grep "_~{chr}.chr10.b38.gmap.gz")
        shapeit4 \
            --input ~{vcf_file} \
            --reference $reference_panel \
            --map $genetic_map \
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
        req_cpu: 16
        #The fixed parameter " req_memory" is used in wdl to specify the running memory
        req_memory: "32Gi"
    }
}
task merge_phased_vcfs_task {
    input {
        #File phase_dir
        Array[File] phased_vcfs
        
    }
    
    #String phased_vcf_files = sep=" " phased_vcfs
    command {
        
        # bcftools concat -o merged_phased.vcf -O v $phased_vcf_files
        bcftools concat -o merged_phased.vcf -O v ${sep=' ' phased_vcfs}
    }
    output {
        File merged_phased_vcf = "merged_phased.vcf"
    }
    runtime {
        docker_url: "public-library/shuliping_b2db8c5c625743c79ca047f52ec29070_public:latest"
        req_cpu: 2
        req_memory: "4Gi"
    }
}

