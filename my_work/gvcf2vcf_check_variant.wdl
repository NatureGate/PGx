version 1.0
workflow PGx{
    input {
        
        #File reference_panel_dir 
        #File genetic_map_dir
        File gvcf_file
        File gvcf_file_index
        String SampleID
        String HLA_dir = "/Files/sz_history/huangfei/BGE/hcWGS/HLA"
        String lush_ref_dir="/Files/sz_history/zhenghaihui/lush"
        String HLA_ref_dir="/Files/sz_history/huangfei/BGE/database"
        String sourceD="megabolt"
        File dbSNPDir = "/Files/sz_history/huangfei/BGE/database/genome/hg38_noalt_withrandom"
        File dbSNP = "/Files/sz_history/huangfei/BGE/database/genome/hg38_noalt_withrandom/All_20180418.new.vcf.gz"
        File dbSNP_index = "/Files/sz_history/huangfei/BGE/database/genome/hg38_noalt_withrandom/All_20180418.new.vcf.gz.csi"
    }
    # String sample_id = basename(vcf_allsite_file, ".allsite.vcf.gz")

    Array[String] chromosomes = ["chr1", "chr2","chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                                 "chr19", "chr20","chr21", "chr22","chrX"]
    
    call gvcf2vcf {
      input:
        infile = gvcf_file,
        infile_index = gvcf_file_index,
        sampleID = SampleID,
        lush_ref_dir = lush_ref_dir

    }
    
    call chrMHandleTask {
        input:
            vcf_allsite_file = gvcf2vcf.vcf,
            vcf_allsite_file_index = gvcf2vcf.vcf_index
    }
    
    scatter (chr in chromosomes) {
        
        call split_vcf{
            input:
                vcf_allsite_file  = gvcf2vcf.vcf,
                vcf_allsite_file_index = gvcf2vcf.vcf_index,
                chr = chr
        }
    }

    call shapeit4_phase_task {
            input:
                vcf_file = split_vcf.vcf_file,
                vcf_file_index = split_vcf.vcf_file_index,
                reference_panel = split_vcf.reference_panel,
                reference_panel_index = split_vcf.reference_panel_index,
                chr_genetic_map = split_vcf.genetic_map,
                threads = 1
    }

    call merge_phased_vcfs_task {
        input:
            
            phased_vcfs = shapeit4_phase_task.phased_vcf,
            vcf_file = split_vcf.vcf_file,
            unphased_vcfs_index = split_vcf.vcf_file_index,
            chrM_vcf = chrMHandleTask.chrM_vcf,
            chrM_vcf_gz = chrMHandleTask.chrM_vcf_gz,
            chrM_vcf_gz_index = chrMHandleTask.chrM_vcf_gz_index,
            dbSNPDir = dbSNPDir
    }

    output{

        File merged_phased_vcf = merge_phased_vcfs_task.merged_phased_vcf
        File merged_unphased_vcf = merge_phased_vcfs_task.merged_unphased_vcf
        File merged_phased_vcf_gz = merge_phased_vcfs_task.merged_phased_vcf_gz
        File merged_phased_vcf_gz_index = merge_phased_vcfs_task.merged_phased_vcf_gz_index
        File merged_unphased_vcf_gz = merge_phased_vcfs_task.merged_unphased_vcf_gz
        File merged_unphased_vcf_gz_index = merge_phased_vcfs_task.merged_unphased_vcf_gz_index

        
    }
}


task gvcf2vcf {
    input {
        String infile
        String infile_index
        String lush_ref_dir
        String sampleID 
    }

    command {
        /usr/local/bin/bcftools annotate -x INFO/MLEAC ${infile} -o ${sampleID}.tmp.vcf
        bgzip ${sampleID}.tmp.vcf && tabix -p vcf ${sampleID}.tmp.vcf.gz
        /home/stereonote/software/gatk-4.2.0.0/gatk --java-options "-Xmx8G -XX:ParallelGCThreads=4" GenotypeGVCFs  -R ${lush_ref_dir}/hg38.fa -O ${sampleID}.allsite.vcf.gz -V ${sampleID}.tmp.vcf.gz -stand-call-conf 10 -all-sites --allow-old-rms-mapping-quality-annotation-data
    }
    runtime {
        req_cpu:32
        req_memory:"64Gi"
        docker_url:"stereonote_hpc/zhenghaihui_5f6bcc1c87c14c54a5a4a82a2bf86fa6_private:latest"


    }
    output {
        File vcf = "${sampleID}.allsite.vcf.gz"
        File vcf_index = "${sampleID}.allsite.vcf.gz.tbi"
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

task chrMHandleTask {
    input {
        File vcf_allsite_file
        File vcf_allsite_file_index
    
        
    }
    

    command <<<
        source activate base
        conda activate bcftools
        bcftools view -o "chrM.allsite.vcf.gz" -O z -r "chrM" "~{vcf_allsite_file}"
        bcftools view -i 'INFO/AC > 0' chrM.allsite.vcf.gz -Oz -o chrM.vcf.gz
        bcftools index chrM.vcf.gz
        bgzip -d -c chrM.vcf.gz > chrM.vcf
    >>>

    output {
        
        File chrM_vcf = "chrM.vcf"
        File chrM_vcf_gz = "chrM.vcf.gz"
        File chrM_vcf_gz_index = "chrM.vcf.gz.csi"
        
    }

    runtime {
        docker_url: "stereonote_hpc/longrui_c6c88994f0244878975c9042eeecb614_private:latest"
        req_cpu: 1
        #The fixed parameter " req_memory" is used in wdl to specify the running memory
        req_memory: "1Gi"
    }
}

task shapeit4_phase_task {
    input {

        Array[File] vcf_file
        Array[File] vcf_file_index
        Array[File] reference_panel
        Array[File] reference_panel_index
        Array[File] chr_genetic_map
        Int threads = 1
        
    }
    # String vcf_file_name = basename(vcf_file, ".vcf.gz")

    command <<<
        echo ~{sep="," vcf_file}
        echo ~{sep="," vcf_file_index}
        echo ~{sep="," chr_genetic_map}
        echo ~{threads}
        echo ~{sep="," reference_panel}
        echo ~{sep="," reference_panel_index}
        all_vcf_file=~{sep="," vcf_file}
        all_vcf_file_index=~{sep="," vcf_file_index}
        all_genetic_map=~{sep="," chr_genetic_map}    
        all_reference_panel=~{sep="," reference_panel}
        all_reference_panel_index=~{sep="," reference_panel_index}
        /opt/conda/bin/python /home/stereonote/script/phase.py --vcf_file ${all_vcf_file}
        
        
    >>>

    output {
        Array[File] phased_vcf = glob("*.phased.vcf")

    }

    runtime {
        docker_url: "stereonote_hpc_external/longrui_8a050cf075d3446ca4a29f06b369c3cc_private:latest"
        req_cpu: 24
        #The fixed parameter " req_memory" is used in wdl to specify the running memory
        req_memory: "32Gi"
    }
}

task merge_phased_vcfs_task {
    input {
        #File phase_dir
        Array[File] phased_vcfs
        Array[File] vcf_file
        Array[File] unphased_vcfs_index
        File chrM_vcf_gz 
        File chrM_vcf_gz_index 
        File chrM_vcf 
        File dbSNPDir
        

    }

    #String phased_vcf_files = sep=" " phased_vcf_files
    command <<<
        source activate base
        conda activate bcftools
        # bcftools concat -o merged_phased.vcf -O v $phased_vcf_files
        bcftools concat -o merged_phased.vcf -O v ~{sep=' ' phased_vcfs} ~{chrM_vcf}
        bcftools concat -o merged_unphased.vcf -O v  ~{sep=' ' vcf_file} ~{chrM_vcf_gz}
        bgzip -c merged_unphased.vcf > unannotate_merged_unphased.vcf.gz
        bgzip -c merged_phased.vcf > unannotate_merged_phased.vcf.gz
        # 添加rsID
        bcftools index unannotate_merged_phased.vcf.gz
        bcftools index unannotate_merged_unphased.vcf.gz
        bcftools annotate -a ~{dbSNPDir}/All_20180418.new.vcf.gz -c ID unannotate_merged_phased.vcf.gz -o merged_phased.vcf.gz
        bcftools annotate -a ~{dbSNPDir}/All_20180418.new.vcf.gz -c ID unannotate_merged_unphased.vcf.gz -o merged_unphased.vcf.gz
        bcftools index merged_phased.vcf.gz
        bcftools index merged_unphased.vcf.gz
    >>>
    output {
        File merged_phased_vcf = "merged_phased.vcf"
        File merged_unphased_vcf = "merged_unphased.vcf"
        File merged_phased_vcf_gz = "merged_phased.vcf.gz"
        File merged_phased_vcf_gz_index = "merged_phased.vcf.gz.csi"
        File merged_unphased_vcf_gz = "merged_unphased.vcf.gz"
        File merged_unphased_vcf_gz_index = "merged_unphased.vcf.gz.csi"
        
    }
    
    runtime {
        docker_url: "stereonote_hpc/longrui_c6c88994f0244878975c9042eeecb614_private:latest"
        req_cpu: 1
        req_memory: "1Gi"
    }
}