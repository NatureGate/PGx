version 1.0
#You need to declaration version information(version 1.0)


workflow PGx{
    input {
        # File vcf_allsite_file
        # File vcf_allsite_file_index
        #File reference_panel_dir 
        #File genetic_map_dir
        File bam_file
        File? hla_file #可能没有输入这个文件
        File bam_index
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
    
     Int len = length(chromosomes)
     Array[Int]  indices = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

    call gvcf2vcf {
      input:
        infile = gvcf_file,
        infile_index = gvcf_file_index,
        sampleID = SampleID,
        lush_ref_dir = lush_ref_dir

    }

    if (!defined(hla_file)) {
        call HLA{
            input:
                sample = SampleID,
                bamFile = select_first([bam_file]),
                bamFile_index = select_first([bam_index]),
                HLA_dir = HLA_dir,
                sourceD = sourceD,
                HLA_ref_dir = HLA_ref_dir

        }
    }

    # 统一拿到最终的 hla.txt
    File finalHla = select_first([hla_file, HLA.hla]) 
   

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
            dbSNPDir = dbSNPDir,
            outfile = cyriusTask.outcall_file,
            sample_id = SampleID
    }
    
    call cyriusTask{
        input:
            bam_file = bam_file,
            hla_file = finalHla,
            bam_index = bam_index
    }

    call preVcfPharmcatTask{
        input:
            vcf_file = merge_phased_vcfs_task.annotated_vcf,
            sample_id = SampleID,
            annotated_vcf = merge_phased_vcfs_task.annotated_vcf
         
    }

    call pharmCatTask{
        input:
            file = preVcfPharmcatTask.pharmcat_vcf,
            outfile = merge_phased_vcfs_task.update_outcall
    }

    call PostPharmCat {
        input:
            report_json_file = pharmCatTask.report_file,
            report_tsv_file = pharmCatTask.report_tsv_result,
            outcall_file = merge_phased_vcfs_task.update_outcall,
            pre_pharmcat_vcf = preVcfPharmcatTask.pharmcat_vcf,
            merged_phased_vcf_gz = merge_phased_vcfs_task.merged_phased_vcf_gz,
            merged_phased_vcf_gz_index = merge_phased_vcfs_task.merged_phased_vcf_gz_index,
            match_json_file = pharmCatTask.match_report,
            sample_id = SampleID,
            merged_unphased_vcf_gz = merge_phased_vcfs_task.merged_unphased_vcf_gz,
            merged_unphased_vcf_gz_index = merge_phased_vcfs_task.merged_unphased_vcf_gz_index
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
        File outcall_file = cyriusTask.outcall_file
        
        File post_pharmcat_file = PostPharmCat.output_file
        File unknown_genotype_file = PostPharmCat.unknown_genotype_file
        File other_file = PostPharmCat.other_file
        File pgx_xls_file = PostPharmCat.pgx_xls_file
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
        bcftools view -o "~{chr}.vcf.gz" -O z -r "~{chr}" "~{vcf_allsite_file}"
        #bcftools view -i 'INFO/AC > 0' ~{chr}.allsite.vcf.gz -Oz -o ~{chr}.vcf.gz
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
        bcftools view -o "chrM.vcf.gz" -O z -r "chrM" "~{vcf_allsite_file}"
        #bcftools view -i 'INFO/AC > 0' chrM.allsite.vcf.gz -Oz -o chrM.vcf.gz
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
        File outfile
        String sample_id
        

    }

    #String phased_vcf_files = sep=" " phased_vcf_files
    command <<<
        source activate base
        conda activate bcftools
        # bcftools concat -o merged_phased.vcf -O v $phased_vcf_files
        bcftools concat -o merged_phased.vcf -O v ~{sep=' ' phased_vcfs} ~{chrM_vcf}
        bcftools concat -o unannotate_merged_unphased.vcf.gz -Oz  ~{sep=' ' vcf_file} ~{chrM_vcf_gz}
        # bgzip -c merged_unphased.vcf > unannotate_merged_unphased.vcf.gz
        bgzip -c merged_phased.vcf > unannotate_merged_phased.vcf.gz
        # 添加rsID
        bcftools index unannotate_merged_phased.vcf.gz
        bcftools index unannotate_merged_unphased.vcf.gz
        bcftools annotate -a ~{dbSNPDir}/All_20180418.new.vcf.gz -c ID unannotate_merged_phased.vcf.gz -o merged_phased.vcf.gz
        bcftools annotate -a ~{dbSNPDir}/All_20180418.new.vcf.gz -c ID unannotate_merged_unphased.vcf.gz -o merged_unphased.vcf.gz
        bcftools index merged_phased.vcf.gz
        bcftools index merged_unphased.vcf.gz
        bcftools filter -e 'FILTER="LowQual"' merged_unphased.vcf.gz -o filtered_merged_unphased.vcf.gz
        bcftools index filtered_merged_unphased.vcf.gz
        bcftools annotate -a merged_phased.vcf.gz -c FORMAT/GT -Ov -o annotated.vcf filtered_merged_unphased.vcf.gz
        bcftools view -i "ID == 'rs267606618' || ID == 'rs267606619' || ID == 'rs267606617' || ID == 'rs28358569'" merged_unphased.vcf.gz|grep -v '^#'|awk -F'\t' '{print $3}'|head -n1 > MTRNR1.txt
        cp ~{outfile}   ~{sample_id}.outcall.tsv
        python /home/stereonote/script/check_MTRNR1.py MTRNR1.txt ~{sample_id}.outcall.tsv
    >>>
    output {
        File merged_phased_vcf = "merged_phased.vcf"
        # File merged_unphased_vcf = "merged_unphased.vcf"
        File merged_phased_vcf_gz = "merged_phased.vcf.gz"
        File merged_phased_vcf_gz_index = "merged_phased.vcf.gz.csi"
        File merged_unphased_vcf_gz = "merged_unphased.vcf.gz"
        File merged_unphased_vcf_gz_index = "merged_unphased.vcf.gz.csi"
        File annotated_vcf = "annotated.vcf"
        File update_outcall = "~{sample_id}.outcall.tsv"
        
        
    }
    
    runtime {
        docker_url: "stereonote_hpc/longrui_16750eaf940e4be288594b3c957f8957_private:latest"
        req_cpu: 1
        req_memory: "2Gi"
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
    java -Xmx1g -jar /home/stereonote/pharmcat-3.0.1-all.jar -vcf ~{file}  -po ~{outfile} -reporterJson --output-dir .  -reporterCallsOnlyTsv

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
        File annotated_vcf
        
      
    }
    
    command <<<
        
        echo ~{vcf_file} ~{annotated_vcf}
        # copy and paste the command to the docker image
        mkdir -p results/pharmcat_ready/
        source activate /home/stereonote/miniconda3
        /home/stereonote/miniconda3/bin/python /home/stereonote/preprocessor/pharmcat_vcf_preprocessor.py \
        -vcf ~{annotated_vcf} \
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
        # File merged_unphased_vcf
        File merged_phased_vcf_gz
        File merged_phased_vcf_gz_index
        File match_json_file
        String sample_id
        File merged_unphased_vcf_gz
        File merged_unphased_vcf_gz_index
    }

    command <<<
        source activate base
        conda activate bcftools
        java -Xmx2g -jar /opt/conda/envs/bcftools/share/snpeff-5.2-1/snpEff.jar -v GRCh38.105  ~{pre_pharmcat_vcf} >annotated.vcf
        bgzip annotated.vcf
        
        bcftools index annotated.vcf.gz
        # bcftools index merged_unphased.vcf.gz
        mkdir -p isec
        bcftools isec  -p isec -Oz -W ~{merged_unphased_vcf_gz}  annotated.vcf.gz
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
         cp ~{sample_id}.pgx.tsv ~{sample_id}.drug.xls
        python /home/stereonote/script/var_drug.py --valid_ids_path valid_ids.txt --sample_id ~{sample_id}
    >>>

    output {
        File output_file = "~{sample_id}.pgx.tsv"
        File unknown_genotype_file = "~{sample_id}_unknown.pgx.tsv"
        File other_file = "~{sample_id}_other.pgx.tsv"
        File pgx_xls_file = "~{sample_id}.drug.xls"
    }

    runtime {
        docker_url: "stereonote_hpc/longrui_0bfae934ea4e484e9897353a4989627c_private:latest"
        req_cpu: 1
        req_memory: "3Gi"
  }
}


task HLA {
    # ancestry analysis
    input{
        # Required parameters
        String sample
        File? bamFile
        File? bamFile_index

        String HLA_dir
        String sourceD
        String HLA_ref_dir
   
    }


    String outdir = sample+"/hla"
    #String out = working_dir + "/" + sample_id + ".int_ancestry_from_vcf.tsv"
    
    String ref = if sourceD == "megabolt" then HLA_ref_dir + "/zboltDB/hg38.fa" else HLA_ref_dir + "/genome/hg38_noalt_withrandom/hg38.fa"
   
    
    command {

        mkdir -p ${outdir}

        if [[ "~{bamFile}" == *.cram ]]; then


            /usr/local/bin/samtools view -b ${bamFile} -T ${ref} -o ${outdir}/${sample}.MHC.bam chr6:28510120-33480577
            /usr/local/bin/samtools view -b -f 4 ${bamFile} -T ${ref} -o ${outdir}/${sample}.unmapped.bam


            /usr/local/bin/samtools view -bh -L ${HLA_dir}/HLA.ref.megabolt.bed -T ${ref} -o ${outdir}/${sample}.HLA.ref.bam ${bamFile}

        else
            /usr/local/bin/samtools view -b ${bamFile}  -o ${outdir}/${sample}.MHC.bam chr6:28510120-33480577
            /usr/local/bin/samtools view -b -f 4 ${bamFile} -o ${outdir}/${sample}.unmapped.bam


            /usr/local/bin/samtools view -bh -L ${HLA_dir}/HLA.ref.megabolt.bed -T ${ref} -o ${outdir}/${sample}.HLA.ref.bam ${bamFile}

        fi

        /usr/local/bin/samtools fastq ${outdir}/${sample}.MHC.bam -0 ${outdir}/${sample}.MHC.fq > ${outdir}/${sample}.MHC.fq
        /usr/local/bin/samtools fastq ${outdir}/${sample}.unmapped.bam -0 ${outdir}/${sample}.unmapped.fq > ${outdir}/${sample}.unmapped.fq


        /usr/local/bin/samtools fastq ${outdir}/${sample}.HLA.ref.bam -0 ${outdir}/${sample}.HLA.ref.fq > ${outdir}/${sample}.HLA.ref.fq


        perl ${HLA_dir}/bin/HLA.pipeline.pl ${outdir} ${outdir}/${sample}

        cp ${outdir}/${sample}/typing.file  ${outdir}/${sample}.hla_typing.txt

    }

    runtime {
        req_cpu: 2
        req_memory: "5Gi"
        docker_url: "stereonote_hpc/zhenghaihui_5f6bcc1c87c14c54a5a4a82a2bf86fa6_private:latest"
    }

    output {
        File hla="${outdir}/${sample}.hla_typing.txt"

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
        req_cpu:4
        req_memory:"8Gi"
        docker_url:"stereonote_hpc/zhenghaihui_5f6bcc1c87c14c54a5a4a82a2bf86fa6_private:latest"


    }
    output {
        File vcf = "${sampleID}.allsite.vcf.gz"
        File vcf_index = "${sampleID}.allsite.vcf.gz.tbi"
    }
}